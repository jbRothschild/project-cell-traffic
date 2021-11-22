import numpy

class Thing():
    def __init__(self, x):
        self.x = x

    def change_x(self, x):
        self.x = x

a = Thing(1)
b = Thing(2)
c = Thing(3)

lst = [a,b,c]
print(lst)

lst2 = [ lst[1], lst[2], lst[0] ]

a.change_x(5)

print(a.x)
print(lst[0].x)
print(lst2[2].x)

lst2[2].change_x(999)

print(a.x)
print(lst[0].x)
print(lst2[2].x)
