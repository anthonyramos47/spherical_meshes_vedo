import vedo 
import numpy as np


# Darboux transformation for a point
def darboux_transform(pt:"np.array[float]", center:"np.array[float]", radius:float, type:int):
    """Sphere inversion of a point pt with respect to a sphere with 
    center and radius.
    """
    # Assert type is 1 or 2
    assert type in [1,2], "Type must be 1 or 2"

    # Distance from pt to the center
    d = np.linalg.norm(np.array(pt) - np.array(center))

    # Sphere inversion
    if d < radius:
        d_P = (radius*radius + (-1)**(type+1) * radius * np.sqrt(radius*radius - d*d)) / d
    elif d > radius:
        d_P = ( - radius*radius + (-1)**(type) * radius * np.sqrt(radius*radius + d*d)) / d
    else: 
        d_P = np.zeros(3)

    return center + d_P* (pt - center) / d

def darboux_transform_mesh_par(mesh:"vedo.Mesh", center:"np.array[float]", radius:float, type:int, clr:str="blue", alph:float=0.7):
    """Sphere inversion of a mesh with respect to a sphere with 
    center and radius.
    """
    # Assert type is 1 or 2
    assert type in [1,2], "Type must be 1 or 2"

    # Get the points and faces
    points = mesh.points()
    faces = mesh.faces()

    # Spherical Triangles
    meshes = []

    # points = [darboux_transform(pt, center, radius, type) for pt in points]

    for f in faces:

        # Get the points of the face
        f_pts = [ points[i] for i in f]

        # Triangulate the face 
        triangles = triangulate_face(f_pts)

        for t in triangles:
            
            pts, faces = spherical_triangle(t[0], t[1], t[2], center, radius)

     
            meshes.append(vedo.Mesh([pts, faces]).color(clr).alpha(alph))

    
    return meshes


def triangulate_face(pts:"np.array[float]"):
    """ Triangulate a polygon with vertices n pts
    """

    if len(pts) != 3:
        bar = np.mean(pts, axis=0)
        
        # Triangles
        triangles = []
        for i in range(len(pts)):
            triangles.append([pts[i], pts[(i+1)%len(pts)], bar])
    else:
        triangles = [pts]

    return triangles


# Unir vector
def unit(v:"np.array[float]"):
    """Unit vector of a vector
    """
    return v / np.linalg.norm(v)

# Function to take the average of all normals of a mesh
def average_normal(mesh:"vedo.Mesh"):
    """Compute the average normal of a mesh
    """
    # Get the points and faces
    points = mesh.points()
    faces = mesh.faces()

    # Compute the normals
    normals = [unit(  np.cross(points[faces[i][1]] - points[faces[i][0]], points[faces[i][2]] - points[faces[i][0]])  )   for i in range(len(faces))]


    # Compute the average normal
    avg_normal = np.mean(normals, axis=0)

    return avg_normal

# Parametriced triangle
def spherical_triangle(p1:"np.array[float]", p2:"np.array[float]", p3:"np.array[float]", center:"np.array[float]", radius:float):
    """Create a parametriced triangle with vertices p1, p2, p3

    return a list of points and a list of faces  
    """
    N = 15

    # U parmeter 
    U = np.linspace(0, 1, N)
    # V parameter
    V = np.linspace(0, 1, N)

    # Create the points

    pts = [ darboux_transform(p1 + u*(p2 - p1) + v*(p3 - p1), center, radius, 1)  \
           for u in U for v in V if u + v <= 1.001]
    
    faces = []

    for u in range(N-1):
        for v in range(N):
            # Auxiliar variable to parametrize the v direction
            nu = sum(range(N, N - v, -1))
            nu1 = sum(range(N, N - (v+1) , -1))

            if u + v < (N-1):
                t1 = [u + nu, (u +1) + nu  , nu1 + u  ] 

                
                faces.append(t1)
                if u +v < N-2:
                    t2 = [(u+1) + nu, (u+1) + nu1, u + (nu1) ]
                
                    faces.append( t2)
            
    return np.array(pts), np.array(faces)




# ---------------------------Quad MESH--------------------------------  
mesh = vedo.load("Models/obj_pq/conical1.obj")

center_DS = np.mean(mesh.points(), axis=0) - 10*average_normal(mesh) 
radius_DS = 25

# Fail triak 
# env1, env2 = darboux_transform_mesh(mesh, center_DS, radius_DS)


DS = vedo.Sphere(center_DS, radius_DS).color('pink').alpha(0.3)

darboux_sph = vedo.Sphere(pos=center_DS, r=radius_DS).color('pink').alpha(0.5)
mesh2 = darboux_transform_mesh_par(mesh, center_DS, radius_DS, 1)

darboux_pts = [ darboux_transform(pt, center_DS, radius_DS, 1) for pt in mesh.points()]
mesh3 = vedo.Mesh([darboux_pts, mesh.faces()]).color('green').alpha(0.8)

merge_mesh = vedo.merge(mesh2).clean()

merge_mesh.subsample(0.00001).wireframe(False)

vedo.write(merge_mesh, "DT_quad_mesh.obj")
vedo.write(mesh, "quad_mesh.obj")
vedo.write(DS, "Darboux_Sphere.obj")

vedo.show(DS, merge_mesh, mesh)

# ---------------------------Test TRIANGLE--------------------------------


# center_DS = np.array([0,0,-10])
# radius_DS = 20

# pts_TM = np.array([[0.1,0,0], [10,10,0], [0,-10,0], [10, 5, -5  ]])
# faces_TM = np.array([[0,2,1], [2,3,1]]) 

# mesh = vedo.Mesh([pts_TM, faces_TM]).color('blue').alpha(0.5)

# enve1, enve2 = darboux_transform_mesh(mesh, center_DS, radius_DS)

# vedo.show(enve1, enve2, mesh, vedo.Sphere(center_DS, radius_DS).color('pink').alpha(0.3))

# ---------------------------TEST QUADS--------------------------------

# center_DS = np.array([0,0,-20])
# radius_DS = 40


# # Points of two quads
# p1 = np.array([-10,0,0])
# p2 = np.array([10,10,0])
# p3 = np.array([0,-10,0])
# p4 = np.array([-10,0,0]) + 1.5*(p2 - p1) + 0.5*(p3 - p1)
# p5 = np.array([10, 17, -2])
# p6 = p4 + 2*(p4 - p2) + 2.5*(p5 - p4)

# pts_TM = np.array([p1, p2, p3, p4, p5, p6])
# faces_TM = np.array([[0,2, 3, 1], [3,5,4,1]])


# mesh = vedo.Mesh([pts_TM, faces_TM]).color('blue').alpha(0.5)

# # label points
# #labels = mesh.labels('id').z(0.1)


# mesh2 = darboux_transform_mesh_par(mesh, center_DS, radius_DS, 1, 'white', 0.8) 

# enve1, enve2 = darboux_transform_mesh(mesh, center_DS, radius_DS)

# DS = vedo.Sphere(center_DS, radius_DS).color('pink').alpha(0.3)

# vedo.show(mesh, mesh2, enve1, DS)





