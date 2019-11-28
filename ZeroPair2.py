class ZeroPair2(object):
    """description of class"""
    def __init__(self, points, color, id):
        self.points = points
        self.color = color
        self.id = id
        self.order = self.__compute_pair_order()
        self.label = self.__compute_label()

        self.used = False

    def __compute_pair_order(self):
        not_null_points = [point for point in self.points if abs(point) != 0]
        return len(not_null_points)
        
    def __compute_label(self):
        label = 'Zero #{0}'.format(self.id)
        return label

