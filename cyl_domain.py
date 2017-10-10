from dolfin import *
import sympy as sp
from sympy.printing import ccode


def expr_body(expr, **kwargs):
    if not hasattr(expr, '__len__'):
        # Defined in terms of some coordinates
        xyz = set(sp.symbols('x[0], x[1], x[2]'))
        xyz_used = xyz & expr.free_symbols
        assert xyz_used <= xyz
        # Expression params which need default values
        params = (expr.free_symbols - xyz_used) & set(kwargs.keys())
        # Body
        expr = ccode(expr).replace('M_PI', 'pi')
        # Default to zero
        kwargs.update(dict((str(p), 0.) for p in params))
        # Convert
        return expr
    # Vectors, Matrices as iterables of expressions
    else:
        return [expr_body(e, **kwargs) for e in expr]


def as_expression(expr, degree=4, **kwargs):
    '''Turns sympy expressions to Dolfin expressions.'''
    return Expression(expr_body(expr), degree=degree, **kwargs)


def mapping_generator(f, bdry_expr, bdry_chi):
    '''TODO'''
    V = f.function_space()
    u = TrialFunction(V)
    v = TestFunction(V)
    
    a = inner(grad(u), grad(v))*dx
    L = inner(Constant(0), v)*dx

    bc = DirichletBC(V, bdry_expr, bdry_chi)

    assembler = SystemAssembler(a, L, bc)
    A = PETScMatrix();
    assembler.assemble(A)

    solver = PETScKrylovSolver('cg', 'hypre_amg')
    solver.set_operators(A, A)
    solver.parameters['relative_tolerance'] = 1E-13
    solver.parameters['absolute_tolerance'] = 1E-13

    x = f.vector()
    b = x.copy()
    while True:
        time_value = yield

        bdry_expr.t = time_value
        assembler.assemble(b)
        
        solver.solve(x, b)

        yield f

        
def compute_mapping(gen, t):
    gen.next()
    return gen.send(t)
        
# -------------------------------------------------------------------

if __name__ == '__main__':
    from numpy.linalg import eigvalsh
    import numpy as np
    import sympy as sp

    A_value = 0.075
    k_value = pi
    c_value = 0.5
    t_value = 0.0

    h_top = Constant(0.4)
    h_bottom = Constant(0.0)
    
    # Variables
    x, y, t = sp.symbols('x[0] x[1] t')
    # Parameters
    A, k, c = sp.symbols('A k c')
    bdry_motion = A*sp.cos(k*y-c*t)**2
    bdry_velocity_x = bdry_motion.diff(t, 1)
    # As expresssions
    bdry_velocity = as_expression((bdry_velocity_x, sp.S(0)),
                                  A=A_value, k=k_value, c=c_value, t=t_value)

    # -------------------------------------------------------------
    
    n = 32
    mesh = RectangleMesh(Point(0.5, -1), Point(1, 1), n, n)

    M = FunctionSpace(mesh, 'CG', 2)
    f = Function(M)

    # Pertuebed Id
    bdry_chi = CompiledSubDomain('near((x[0]-0.5)*(x[0]-1), 0)')
    # bdry_chi = DomainBoundary()
    bdry_expr = Expression('near(x[0], 0.5) ? %s : 0' % ccode(bdry_motion).replace('M_PI', 'pi'),
                           A=A_value, k=k_value, c=c_value, t=t_value, degree=3)

    top, bottom = 2, 1 
    top_chi = CompiledSubDomain('near(x[1], 1)')
    bottom_chi = CompiledSubDomain('near(x[1], -1)')

    bdries = FacetFunction('size_t', mesh, 0)
    top_chi.mark(bdries, top)
    bottom_chi.mark(bdries, bottom)
    

    # ----------------------------------------------------------------
    
    v_elm = VectorElement('Lagrange', triangle, 2)
    p_elm = FiniteElement('Lagrange', triangle, 1)
    TH = MixedElement([v_elm, p_elm])
    W = FunctionSpace(mesh, TH)

    u, p = TrialFunctions(W)
    v, q = TestFunctions(W)

    x, y = SpatialCoordinate(mesh)
    n = FacetNormal(mesh)
    # The mapping
    r = x+f
    fmap = as_vector((r, y))
    F = grad(fmap)
    J = det(F)

    Grad = lambda arg: dot(grad(arg), inv(F))
    Div = lambda arg: inner(grad(arg), inv(F))

    a = inner(Grad(u), Grad(v))*J*r*dx + inner(u[0]/r, v[0]/r)*J*r*dx +\
        inner(p, Div(v))*J*r*dx + inner(p, v[0]/r)*J*r*dx +\
        inner(q, Div(u))*J*r*dx + inner(q, u[0]/r)*J*r*dx
    
    L = inner(Constant((0, 0)), v)*J*r*dx + \
        h_top*inner(v, dot(transpose(inv(F)), n))*J*ds(top) +\
        h_bottom*inner(v, dot(transpose(inv(F)), n))*J*ds(bottom)
              
    bcs = [DirichletBC(W.sub(0), bdry_velocity, 'near(x[0], 0.5)'),
           DirichletBC(W.sub(0), Constant((0, 0)), 'near(x[0], 1)')]

    wh = Function(W)

    V = FunctionSpace(mesh, v_elm)
    Q = FunctionSpace(mesh, p_elm)
    
    uh, ph = Function(V), Function(Q)
    
    assigner = FunctionAssigner([V, Q], W)
    
    # ---------------------------------------------------------------
              
    xy = mesh.coordinates().copy()
    x, y = xy.T

    f_mesh = Mesh(mesh)
    fuh = Function(FunctionSpace(f_mesh, v_elm))
    fph = Function(FunctionSpace(f_mesh, p_elm))
    
    gen = mapping_generator(f, bdry_expr, bdry_chi)

    t_min = []
    u_out, p_out = File('./data/F_uh.pvd'), File('./data/F_ph.pvd')

    S = FunctionSpace(mesh, 'DG', 0)
    for t in np.linspace(0, 5, 200):
        f = compute_mapping(gen, t)
        bdry_velocity.t = t
       
        A, b = assemble_system(a, L, bcs)

        # print np.sort(np.abs(eigvalsh(A.array())))
        
        solve(A, wh.vector(), b)

        assigner.assign([uh, ph], wh)

        # Modif f mesh
        new = np.array([f(xi, yi) for xi, yi in zip(x, y)])
        f_mesh.coordinates()[:, 0] = x + new

        fuh.vector().zero(); fuh.vector().axpy(1, uh.vector())
        fph.vector().zero(); fph.vector().axpy(1, ph.vector())
        

        u_out << (fuh, t)
        p_out << (fph, t)

        print t, wh.vector().norm('l2')
        t_min.append(min(cell.volume() for cell in cells(f_mesh)))
    print 'Tmin', min(t_min)
