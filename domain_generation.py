import subprocess, os, shutil
from dolfin import Mesh, MeshFunction


def gmsh_closed_polygon(points, line_tags=None, size=1.):
    '''
    Closed polygon defined by seguence of points (x, y, s) where
    s*size will be the characteristic size. Line tags is an array giving
    each line on the boundary a value.
    '''
    geo = 'domain.geo'
    with open(geo, 'w') as f:
        f.write('SIZE = %g;\n' % size)
        
        point = 'Point(%d) = {%.16f, %.16f, 0, %g*SIZE};\n'
        for pindex, p in enumerate(points, 1):
            f.write(point % ((pindex, ) + p))

        line = 'Line(%d) = {%d, %d};\n'
        for lindex in range(1, len(points)):
            f.write(line % (lindex, lindex, lindex+1))
        # Close loop
        nlines = lindex + 1
        f.write(line % (nlines, nlines, 1))

        # Isolate unique tags
        if line_tags is None: line_tags = [1]*nlines

        unique_tags = set(line_tags)
        for tag in unique_tags:
            indices = filter(lambda i: line_tags[i] == tag, range(nlines))
            indices = ','.join(['%d' % (index + 1) for index in indices])
            
            f.write('Physical Line(%d) = {%s};\n' % (tag, indices))

        loop = ','.join(['%d' % l for l in range(1, nlines+1)])
        f.write('Line Loop(1) = {%s};\n' % loop)
        f.write('Plane Surface(1) = {1};\n')
        f.write('Physical Surface(1) = {1};\n')
    return geo


def geo_to_xml(geo, scale, save='', dir='./meshes'):
    '''Gmsh -> dolfin-convert (optinally save)'''
    root, _ = os.path.splitext(geo)
    msh_file = '.'.join([root, 'msh'])
    subprocess.call(['gmsh -2 -clscale %g -optimize %s' % (scale, geo)], shell=True)
    assert os.path.exists(msh_file)

    # Convert to xdmf
    xml_file = '%s.xml'
    xml_facets = '%s_facet_region.xml'
    xml_volumes = '%s_physical_region.xml'

    subprocess.call(['dolfin-convert %s %s' % (msh_file, xml_file % root)], shell=True)
    # All 3 xml files should exist
    assert all(os.path.exists(f % root) for f in (xml_file, xml_facets, xml_volumes))

    if save:
        assert os.path.exists(dir) and os.path.isdir(dir)

        for f in (xml_file, xml_facets):
            shutil.move(f % root, os.path.join(dir, f % save))

        root = save
    else:
        dir = '.'
        
    return [os.path.join(dir, xml_file % root), os.path.join(dir, xml_facets % root)]


def compute_polygon(points):
    ''''
    You are given p0-------p5
                   \        \p4
                    \      /
                    p1    |
                    /     |
                  p2-------p3

    Here I insert point to left and right side.
    '''
    assert len(points) == 6
    # For same height the symmetry is already there
    if abs(points[1][1] - points[4][1]) < 1E-13:
        # Top is 1, bottom is 2. Left takes first
        tags = [3, 4, 2, 5, 6, 1]
        return points, tags

    def insert_point(p, (A, B)):
        s = (A[1]-p[1])/float(A[1]-B[1])
        #
        return (A[0] + s*(B[0]-A[0]), A[1] + s*(B[1]-A[1]), p[-1])
    
    # If p1 is higher it inserts into p5, p4 and p4 is inserted p1 and p2
    if points[1][1] > points[4][1]:
        points = (points[0],
                  points[1],
                  insert_point(points[4], (points[1], points[2])),
                  points[2],
                  points[3],
                  points[4],
                  insert_point(points[1], (points[4], points[5])),
                  points[5])
    else:
        points = (points[0],
                  insert_point(points[4], (points[0], points[1])),
                  points[1],
                  points[2],
                  points[3],
                  insert_point(points[1], (points[3], points[4])),
                  points[4],
                  points[5])
        
    tags = [3, 4, 5, 2, 6, 7, 8, 1]
    return points, tags


def generate_mesh(r_inner, r_outer, length, inner_p=None, outer_p=None, inner_size=1., outer_size=1.,
                  size=1., scale=0.5, save='', dir='./meshes'):
    '''Special case for peristalsis'''
    # Try loading
    if save:
        save = '_'.join([save, str(scale)])
        xml_file = os.path.join(dir, '%s.xml' % save)
        xml_facets = os.path.join(dir, '%s_facet_region.xml' % save)

        xmls = [xml_file, xml_facets]
        if all(os.path.exists(xml) for xml in xmls):
                mesh = Mesh(xmls[0])
                boundaries = MeshFunction('size_t', mesh, xmls[1])

                return mesh, boundaries

    # Generate
    assert r_outer > r_inner > 0
    assert length > 0
    # Setup geo

    # Left half
    points = [(r_inner, length/2., inner_size)]
    if inner_p is not None:
        assert -length/2. < inner_p[1] < length/2.
    else:
        inner_p = (r_inner, 0)
    points.append(inner_p + (inner_size, ))
        
    points.append((r_inner, -length/2., inner_size))

    # Right half
    points.append((r_outer, -length/2., outer_size))
    
    if outer_p is not None:
        assert -length/2. < outer_p[1] < length/2.
    else:
        outer_p = (r_outer, 0)
    points.append(outer_p + (outer_size, ))

    points.append((r_outer, length/2., outer_size))

    # Fill in the extrapoints 
    polygon, tags = compute_polygon(points)
    geo = gmsh_closed_polygon(polygon, line_tags=tags, size=1.)
    xmls = geo_to_xml(geo, scale, save=save)

    # Return as the output that peristalsis solver will use
    mesh = Mesh(xmls[0])
    boundaries = MeshFunction('size_t', mesh, xmls[1])

    return mesh, boundaries

# -------------------------------------------------------------------

if __name__ == '__main__':
    from dolfin import plot, interactive
    
    mesh, bdries = generate_mesh(r_inner=0.5,
                                 r_outer=1.0,
                                 length=1,
                                 inner_p=(0.6, -0.2),
                                 outer_p=(0.9, 0.),
                                 inner_size=0.5,
                                 outer_size=1.,
                                 size=1.,
                                 scale=1./2**5,
                                 save='')

    plot(bdries)
    interactive()
    
