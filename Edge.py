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




# Example Usage:

edge1 = Edge([0, 0, 0], [1, 1, 1])
edge2 = Edge((0, 0, 0), (2, 2, 2))
edge3 = Edge((0, 0, 0), (1, 1, 1))  # This will be considered equal to edge1

