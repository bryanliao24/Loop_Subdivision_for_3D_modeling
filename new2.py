import mesh
import math
import numpy as np
from collections import defaultdict

class Edge:
    def __init__(self, vertex1, vertex2):
        self.vertices = (vertex1, vertex2)
        self.vertices = tuple(sorted(self.vertices))

    def __eq__(self, other):
        return self.vertices == other.vertices

    def __hash__(self):
        return hash(self.vertices)

    def __repr__(self):
        return f"Edge({self.vertices[0]}, {self.vertices[1]})"

def is_boundary(halfedge):
    if halfedge.opposite == None:
        return True
    else:
        return False

def create_edge(vertex_index1, vertex_index2):
    vertex1 = vertices[vertex_index1].get_vertex()
    vertex2 = vertices[vertex_index2].get_vertex()
    edge = Edge(vertex1, vertex2)
    return edge

def get_coordinate(vertex_index):
    return vertices[vertex_index].get_vertex()

def sharedFace(my_edge, face_points):
    remote_points = [0] * 2
    for fp in face_points:
        if my_edge.vertices[0] != fp and my_edge.vertices[1] != fp:
            remote_points = fp

    return remote_points

def isShareEdge(my_edge, facets):
    count = 0
    remote_point_list = []
    for face in facets:
        a = get_coordinate(face.a)
        b = get_coordinate(face.b)
        c = get_coordinate(face.c)

        face_points = [a,b,c]
        
        if my_edge.vertices[0] in face_points and my_edge.vertices[1] in face_points:
            remote_point_list.append(sharedFace(my_edge, face_points))
            count += 1
        if count == 2:
            # print("shared edge", count)
            return remote_point_list
    if count == 1:
        # print("boundary edge",count)
        # this is boundary
        return remote_point_list
    else:
        return

# update all the old interior vertex
def get_interior_odd_vertex(edge, remote1, remote2):
    v0 = edge.vertices[0]
    v1 = remote1
    v2 = edge.vertices[1]
    v3 = remote2

    new_x = 3/8 * (v0[0] + v2[0]) + 1/8 * (v3[0] + v1[0])
    new_y = 3/8 * (v0[1] + v2[1]) + 1/8 * (v3[1] + v1[1])
    new_z = 3/8 * (v0[2] + v2[2]) + 1/8 * (v3[2] + v1[2])
    return [new_x, new_y, new_z]

# update all the old boundary edge vertex
def get_boundary_odd_vertex(edge):
    v0 = edge.vertices[0]
    v2 = edge.vertices[1]

    new_x = 1/2 * (v0[0] + v2[0])
    new_y = 1/2 * (v0[1] + v2[1])
    new_z = 1/2 * (v0[2] + v2[2])
    
    return [new_x, new_y, new_z]


bianjie = defaultdict(list)
boundary = []
def get_odd_vertice_and_new_facets(vertices, facets):
    new_facets = []
    new_vertices = []
    
    for vertex in vertices:
        new_vertices.append(vertex.get_vertex())

    for face in facets:

        a = get_coordinate(face.a)
        b = get_coordinate(face.b)
        c = get_coordinate(face.c)

        edge1 = Edge(a,b)
        edge2 = Edge(b,c)
        edge3 = Edge(c,a)

        remote_points_edge_1 = isShareEdge(edge1, facets)
        if len(remote_points_edge_1) == 2:
            odd_vertex_ab = get_interior_odd_vertex(edge1, remote_points_edge_1[0], remote_points_edge_1[1])
        elif len(remote_points_edge_1) == 1: 
            # boundary mapping
            odd_vertex_ab = get_boundary_odd_vertex(edge1)
            boundary.append(face.a)
            bianjie[face.a].append(face.b)
            bianjie[face.b].append(face.a)
    
        remote_points_edge_2 = isShareEdge(edge2, facets)
        if len(remote_points_edge_2) == 2:
            odd_vertex_bc = get_interior_odd_vertex(edge2, remote_points_edge_2[0], remote_points_edge_2[1])
        elif len(remote_points_edge_2) == 1: 
            odd_vertex_bc = get_boundary_odd_vertex(edge2)
            boundary.append(face.b)
            bianjie[face.c].append(face.b)
            bianjie[face.b].append(face.c)

        remote_points_edge_3 = isShareEdge(edge3, facets)
        if len(remote_points_edge_3) == 2:
            odd_vertex_ca = get_interior_odd_vertex(edge3, remote_points_edge_3[0], remote_points_edge_3[1])
        elif len(remote_points_edge_3) == 1:  
            odd_vertex_ca = get_boundary_odd_vertex(edge3)
            boundary.append(face.c)
            bianjie[face.c].append(face.a)
            bianjie[face.a].append(face.c)
 
        if odd_vertex_ab in new_vertices:
            index_ab = new_vertices.index(odd_vertex_ab)
        else:
            new_vertices.append(odd_vertex_ab)
            index_ab = new_vertices.index(odd_vertex_ab)

        if odd_vertex_bc in new_vertices:
            index_bc = new_vertices.index(odd_vertex_bc)
        else:
            new_vertices.append(odd_vertex_bc)
            index_bc = new_vertices.index(odd_vertex_bc)

        if odd_vertex_ca in new_vertices:
            index_ca = new_vertices.index(odd_vertex_ca)
        else:
            new_vertices.append(odd_vertex_ca)
            index_ca = new_vertices.index(odd_vertex_ca)     
        
        new_face_0 = [face.a, index_ab, index_ca]
        new_face_1 = [index_ab, face.b, index_bc]
        new_face_2 = [index_bc, face.c, index_ca]
        new_face_3 = [index_ab, index_bc, index_ca]

        new_facets.append(new_face_0)
        new_facets.append(new_face_1)
        new_facets.append(new_face_2)
        new_facets.append(new_face_3)

    return new_vertices, new_facets

def calculate_neighbor(vertices, facets):
    neighbor_index = defaultdict(list) # this is all the mapping with one vertex to all neighbor index 
    neighbor_num = {} # this is all the mapping with one vertex to its number of neighbors 

    for vertex in vertices:

        for f in facets:
            fff = [f.a, f.b, f.c]
            if vertex.index in fff:
                neighbor_num[vertex.index] = neighbor_num.get(vertex.index, 0) + 1
                for i in fff:
                    if vertex.index != i and i not in neighbor_index[vertex.index]:
                        neighbor_index[vertex.index].append(i)
    # print(neighbor_index)       
    return neighbor_index, neighbor_num

def calculate_beta(neighbors_count):
    if neighbors_count == 3:
        return 3/16
    else:
        return (1.0 / neighbors_count) * (5.0 / 8.0 - ((3.0 / 8.0) + (1.0 / 4.0) * math.cos(2 * math.pi / neighbors_count))** 2)
# print(len(neighbor_num)) -> 8 

def Kclosest(arr, x, k):
    # if x lies beyond arr[0]
    if x < arr[0]:
        arr = arr[:k]
        return arr
    # if x lies beyond arr[-1]
    elif x > arr[-1]:
        arr = arr[len(arr)-k:]
        return arr
    # Apply sliding window on the aux array
    elif x >= arr[0] and x <= arr[-1]:
        aux = []
        for ele in arr:
            diff = int(math.fabs(ele - x))
            aux.append(diff)
        i = 0
        j = 0
        s = 0
        maxi_s = float("inf")
        # find window with smallest sum
        while j < len(aux):
            s += aux[j]
            if (j - i + 1) < k:
                j += 1
            elif (j - i + 1) == k:
                if s < maxi_s:
                    maxi_s = s
                    start = i
                    end = j
                s -= aux[i]
                i += 1
                j += 1
        # Just reutrn the first window with smallest sum 
        arr = arr[start:end+1]
        return arr
    
import collections
dictt = collections.defaultdict(list) 
def update_old_vertice(vertices, facets):
    neighbor_index, neighbor_num = calculate_neighbor(vertices, facets)

    bound = [] # store all the boundary index into list
    count_interval = 0
    count_allsmall = 0
    count_allbig = 0

    bbo = set(boundary)
   
    update_old = []
    for i in range(len(neighbor_num)):

        total_vi_x = 0
        total_vi_y = 0
        total_vi_z = 0

        if i in bbo:
            # print("yes")
            # arr = neighbor_index[i]
            # arr.sort()
            # # 初始化最接近i的两个数字
            # # close2 = Kclosest(arr, i, 2)
            # # closest_num1 = close2[0]
            # # closest_num2 = close2[1]
            # # 找到最大的数
            # if i < 240: 
            #     max_number = max(arr)

            #     # 删除最大的数
            #     arr.remove(max_number)
            # closest_num1 = None
            # closest_num2 = None
            # # 遍历数组找到第一个大于等于i的数字
            # for num in arr:
            #     if num >= i:
            #         closest_num1 = num
            #         break

            # # 如果找到了第一个大于等于1的数字，找到它的前一个数字
            # if closest_num1 is not None:
            #     index = arr.index(closest_num1)
            #     if index > 0:
            #         closest_num2 = arr[index - 1]
            #         count_interval += 1
            # if closest_num1 is None or closest_num2 is None:
            #     close2 = Kclosest(arr, i, 2)
            #     closest_num1 = close2[0]
            #     closest_num2 = close2[1]
            # if closest_num1 is None or closest_num2 is None:# 全大全小情况
            #     target = i
            #     # 初始化标志变量
            #     all_greater = True
            #     all_smaller = True
                
            #     # 遍历数组并检查每个元素
            #     for num in arr:
            #         if num <= target:
            #             # 如果遇到一个元素不大于给定数字，将 all_greater 设置为False
            #             all_greater = False
            #         if num >= target:
            #             # 如果遇到一个元素不小于给定数字，将 all_smaller 设置为False
            #             all_smaller = False
            #     if all_greater:
            #         closest_num1 = arr[0]
            #         closest_num2 = arr[1]
            #         count_allbig += 1
            #     if all_smaller:
            #         closest_num1 = arr[-1]
            #         closest_num2 = arr[-2]
            #         count_allsmall += 1
            closest_num1 = bianjie[i][0]
            closest_num2 = bianjie[i][1]
            dictt[i].append(closest_num1)
            dictt[i].append(closest_num2)
            # print(i)
            total_vi_x = vertices[closest_num1].get_vertex()[0] + vertices[closest_num2].get_vertex()[0] 
            total_vi_y = vertices[closest_num1].get_vertex()[1] + vertices[closest_num2].get_vertex()[1]
            total_vi_z = vertices[closest_num1].get_vertex()[2] + vertices[closest_num2].get_vertex()[2]

            v_update_x = (3.0 / 4.0) * vertices[i].get_vertex()[0] + (1.0 / 8.0) * total_vi_x
            v_update_y = (3.0 / 4.0) * vertices[i].get_vertex()[1] + (1.0 / 8.0) * total_vi_y
            v_update_z = (3.0 / 4.0) * vertices[i].get_vertex()[2] + (1.0 / 8.0) * total_vi_z
            update_old.append([round(v_update_x, 7), round(v_update_y, 7), round(v_update_z, 7)])
            #update_old.append([v_update_x, v_update_y, v_update_z])
        else:
            for j in neighbor_index[i]:
                total_vi_x += calculate_beta(neighbor_num[i]) * vertices[j].get_vertex()[0] 
            for j in neighbor_index[i]:
                total_vi_y += calculate_beta(neighbor_num[i]) * vertices[j].get_vertex()[1] 
            for j in neighbor_index[i]:
                total_vi_z += calculate_beta(neighbor_num[i]) * vertices[j].get_vertex()[2] 
            v_update_x = (1-neighbor_num[i]*calculate_beta(neighbor_num[i])) * vertices[i].get_vertex()[0] + total_vi_x
            v_update_y = (1-neighbor_num[i]*calculate_beta(neighbor_num[i])) * vertices[i].get_vertex()[1] + total_vi_y
            v_update_z = (1-neighbor_num[i]*calculate_beta(neighbor_num[i])) * vertices[i].get_vertex()[2] + total_vi_z

            update_old.append([round(v_update_x, 6), round(v_update_y, 6), round(v_update_z, 6)])
            #update_old.append([v_update_x, v_update_y, v_update_z])
    # print(neighbor_num)
    # print(neighbor_index)
    # print(count_interval)
    # print(count_allsmall)
    # print(count_allbig)

    return update_old


def save_halfmesh_as_obj(all_vertex, new_facet, file_name):
    with open(file_name, 'w') as open_file:
        for vertex in all_vertex:
            x = vertex[0]
            y = vertex[1]
            z = vertex[2]
            open_file.write("v {} {} {}\n".format(x, y, z))

        for face in new_facet:
            f0 = face[0]
            f1 = face[1]
            f2 = face[2]
            open_file.write("f {} {} {}\n".format(f0+1, f1+1, f2+1))

data_path = 'tests/data/bunny.off'

# HalfedgeMesh
mesh = mesh.HalfedgeMesh(data_path)

vertices, facets = mesh.vertices, mesh.facets
new_vertice, new_facets = get_odd_vertice_and_new_facets(vertices, facets)

old_updated_vertice = update_old_vertice(vertices, facets)

new_vertice[:(len(old_updated_vertice))] = old_updated_vertice

save_halfmesh_as_obj(new_vertice, new_facets, 'iteration_1.obj')


print(boundary)
# print(dictt) # 程式实际选的边界邻居点
# print(len(bianjie))
