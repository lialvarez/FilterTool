import numpy as np

class PolePair2(object):
    """description of class"""
    def __init__(self, points, color, id):
        self.points = points
        self.color = color
        self.id = id
        self.order = self.__compute_pair_order()
        self.Q = self.__compute_q()
        self.label = self.__compute_label()

        self.used = False

    def __compute_pair_order(self):
        not_null_points = [point for point in self.points if abs(point) != 0]
        return len(not_null_points)
    
    def __compute_q(self):
        # TODO: only case considered conj pair:
        re = np.real(self.points[0])
        mag = abs(self.points[0])
        return abs((mag/(2*re)))

    def __compute_label(self):
        label = 'Pole #{0} - Q: {1:.2f}'.format(self.id, self.Q)
        return label

