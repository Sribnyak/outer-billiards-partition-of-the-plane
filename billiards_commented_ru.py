from fractions import Fraction as Frac  # поддержка рациональных чисел
import matplotlib.pyplot as plt  # библиотека для построения графиков


TOO_MANY_BORDERS = 70000


# Классы для хранения векторов, отрезков и пр.

class FracVec2:
    def __init__(self, x, y):
        self.x = Frac(x)
        self.y = Frac(y)

    def __add__(self, other):
        return FracVec2(self.x + other.x, self.y + other.y)
    def __sub__(self, other):
        return FracVec2(self.x - other.x, self.y - other.y)
    def __neg__(self):
        return FracVec2(-self.x, -self.y)
    def __rmul__(self, num):
        return FracVec2(self.x * num, self.y * num)
    def crs(self, other):
        return self.x * other.y - self.y * other.x
    def dot(self, other):
        return self.x * other.x + self.y * other.y
    def __eq__(self, other):
        return (type(other) == type(self) and
                self.x == other.x and self.y == other.y)

    def reflect(self, other):  # отражает объект other относительно точки self
        if isinstance(other, Segment):
            return Segment(self.reflect(other.begin), self.reflect(other.end))
        if isinstance(other, Ray):
            return Ray(self.reflect(other.begin), -other.vector)
        return self - (other - self)


class Segment:
    def __init__(self, begin, end):
        self.begin = begin
        self.end = end

    def __contains__(self, point):
        return (self.begin - point).dot(self.end - point) < 0

    def as_ray(self):  # необходимо для нахождения пересечений
        return Ray(self.begin, self.end - self.begin)

    def split_xy(self):  # необходимо для построения картинки в matplotlib
        return [self.begin.x, self.end.x], [self.begin.y, self.end.y]

    def max_coord(self):
        return max(max(abs(p.x), abs(p.y)) for p in [self.begin, self.end])


class Ray:
    def __init__(self, begin, vector):
        self.begin = begin
        self.vector = vector

    def __contains__(self, point):
        vector = point - self.begin
        return vector.crs(self.vector) == 0 and vector.dot(self.vector) >= 0

    def intersection(self, other):
        if self.vector.crs(other.vector) == 0:
            return None
        t = ((other.begin - self.begin).crs(other.vector) /
             self.vector.crs(other.vector))
        p = self.begin + t * self.vector
        if t > 0 and p in other:
            return p
        return None

    def as_ray(self):  # необходимо для нахождения пересечений
        return self

    def cut(self, max_coord):  # превращает в отрезок - для построения картинки
        if self.vector.x == 0:
            self.cut_y(max_coord)
        else:
            self.cut_x(max_coord)
            if (self.begin.y + self.vector.y > max_coord or
                    self.begin.y + self.vector.y < -max_coord):
                self.cut_y(max_coord)
        return Segment(self.begin, self.begin + self.vector)

    def cut_x(self, max_coord):
        if self.vector.x > 0:
            k = (max_coord - self.begin.x) / self.vector.x
        else:
            k = -(max_coord + self.begin.x) / self.vector.x
        self.vector = k * self.vector

    def cut_y(self, max_coord):
        if self.vector.y > 0:
            k = (max_coord - self.begin.y) / self.vector.y
        else:
            k = -(max_coord + self.begin.y) / self.vector.y
        self.vector = k * self.vector

    def max_coord(self):
        return max(max(abs(p.x), abs(p.y))
                   for p in [self.begin, self.begin + self.vector])


class Zone:  # для хранения областей, для которых T^-1 действует как движение
    def __init__(self, corner, left_vector, right_vector):
        self.corner = corner
        self.left_vector = left_vector
        self.right_vector = right_vector
        self.prv = None  # для хранения соседних областей
        self.nxt = None

    def __contains__(self, point):
        v = point - self.corner
        return v.crs(self.left_vector) > 0 and v.crs(self.right_vector) < 0

    def get_left_ray(self):
        return Ray(self.corner, self.left_vector)

    def get_right_ray(self):
        return Ray(self.corner, self.right_vector)

    def intersection(self, ray):
        p = ray.intersection(self.get_left_ray())
        if p is not None:
            return p
        return ray.intersection(self.get_right_ray())


class Polygon:
    def __init__(self, vertices):
        self.size = len(vertices)
        self.vertices = vertices

    def get_edges(self):
        return [Segment(self.vertices[i-1], self.vertices[i])
                for i in range(self.size)]


# Главный класс с основной логикой

class Field:
    def __init__(self, table):
        self.table = table
        self.zones = []
        for i in range(table.size):
            corner = table.vertices[i-1]
            left = table.vertices[i-2] - corner
            right = corner - table.vertices[i]
            self.zones.append(Zone(corner, left, right))
        for i in range(table.size):
            self.zones[i].prv = self.zones[i-1]
            self.zones[i-1].nxt = self.zones[i]
        self.iteration = 0
        self.borders = []

    def find_borders(self, max_iteration):
                       # храним пары (граница, итерация)
        self.borders = [(Ray(edge.end, edge.end - edge.begin), 1)
                        for edge in self.table.get_edges()]
        # итерация здесь - первая не определённая итерация для точек границы

        used_borders = 0  # храним, для скольки границ вычислена T^-1

        while used_borders < len(self.borders):
            border, iteration = self.borders[used_borders]
            if iteration != self.iteration:  # переходим к новой итерации
                self.iteration = iteration
                if iteration % 5 == 0:
                    print(f'Iteration {iteration}')
                if iteration == max_iteration:
                    print(f'Finishing on iteration {iteration}')
                    print()
                    break
            self.borders += [(zone.corner.reflect(border_part), iteration + 1)
                             for border_part, zone in self.split_border(border)]
            used_borders += 1
            if len(self.borders) == TOO_MANY_BORDERS:
                max_iteration = iteration + 1
                print(f'Got too many borders.')
    
    def split_border(self, border, begin_zone=None):
        if begin_zone is None:  # при рекурсивном вызове begin_zone известна
            begin_zone = self.get_begin_zone(border)
        if begin_zone is None:
            return []
        intersection = begin_zone.intersection(border.as_ray())
        if intersection is None:
            return [(border, begin_zone)]
        if isinstance(border, Ray):
            border1 = Segment(border.begin, intersection)
            border2 = Ray(intersection, border.vector)
        elif intersection in border:
            border1 = Segment(border.begin, intersection)
            border2 = Segment(intersection, border.end)
        else:
            return [(border, begin_zone)]
        if intersection in begin_zone.get_left_ray():
            new_zone = begin_zone.prv
        else:
            new_zone = begin_zone.nxt
        return [(border1, begin_zone)] + self.split_border(border2, new_zone)
    
    def get_begin_zone(self, border):
        for zone in self.zones:
            if border.begin in zone:
                return zone
            if (border.begin in zone.get_left_ray() and
                    border.begin != zone.corner):
                product = border.as_ray().vector.crs(zone.left_vector)
                if product < 0:
                    return zone.prv
                elif product > 0:
                    return zone
        return None

    def show(self, colour='black', linewidth=0.5):
        field_size = max(border.max_coord() for border, _ in self.borders) * 2
        x = [p.x for p in self.table.vertices]
        y = [p.y for p in self.table.vertices]
        plt.plot(x + [x[0]], y + [y[0]], colour, linewidth=linewidth)
        plt.fill(x, y, colour)
        n = len(self.borders)
        for i, (border, _) in enumerate(self.borders):
            if isinstance(border, Ray):
                x, y = border.cut(field_size).split_xy()
            else:
                x, y = border.split_xy()
            plt.plot(x, y, colour, linewidth=linewidth)
            progress = i * 100 // n
            if progress % 5 == 0 and (i - 1) * 100 // n < progress:
                print(f'{progress}% of the image is ready')
        plt.title('Set $B_{'+str(self.iteration)+'}$')
        plt.show()


def main():
    print('Enter coordinates of your polygon clockwise in the format x1 y1 x2 y2 ...')
    print('For example: 0 0 -1/2 1 1/2 3/2 3/2 1 1 0')
    print('Another example: 3 1 3 -1 1 -3 -1 -3 -3 -1 -3 1 -1 3 1 3')
    print('Or if your polygon is a trapezium, enter only its base ratio')
    print('For example: 4/7')

    polygon = None
    while polygon is None:
        inp = input('Your input: ').split()
        if len(inp) == 1:  # трапеция
            polygon = Polygon([FracVec2(1, 0), FracVec2(0, 0),
                               FracVec2(0, 1), FracVec2(Frac(inp[0]), 1)])
        elif inp and len(inp) % 2 == 0:
            polygon = Polygon([FracVec2(Frac(inp[i*2]), Frac(inp[i*2+1]))
                               for i in range(len(inp)//2)])

    max_iteration = int(input('Max iteration: '))
    field = Field(polygon)
    field.find_borders(max_iteration)
    field.show()

main()
