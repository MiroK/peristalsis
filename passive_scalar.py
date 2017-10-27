from dolfin import *

mesh = RectangleMesh(Point(-1, -1), Point(1, 1), 16, 16)

V = FunctionSpace(mesh, 'DG', 0)
u = TrialFunction(V)
v = TestFunction(V)

dt = Constant(1E-3)

bdries = FacetFunction('size_t', mesh, 0)
CompiledSubDomain('near(x[0], 0)').mark(bdries, 1)

ds = Measure('ds', subdomain_data=bdries)

w = Expression(('-x[1]', 'x[0]'), degree=1)
u0 = interpolate(Expression('(pow(x[0]-0.25, 2) + pow(x[1]-0.25, 2)) < pow(0.25, 2) ? 1 : 0', degree=2), V)

n = FacetNormal(mesh)

flux = lambda u, w, n: dot(avg(w), n('+'))*avg(u) + Constant(0.5)*abs(dot(avg(w), n('+')))*jump(u)

a = inner(u, v)*dx + Constant(dt/2)*inner(flux(u, w, n), jump(v))*dS
L = inner(u0, v)*dx - Constant(dt/2)*inner(flux(u0, w, n), jump(v))*dS

A = assemble(a)
b = assemble(L)

out = File('scalar.pvd')

T = 0
while T < 1:
    T += dt(0)

    assemble(L, b)
    solve(A, u0.vector(), b)

    out << u0, T
