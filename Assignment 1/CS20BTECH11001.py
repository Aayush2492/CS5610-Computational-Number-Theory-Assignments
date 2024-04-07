import csv

data = []

with open('input-gcd.csv', 'r') as f:
    reader = csv.reader(f, delimiter=',')
    for row in reader:
        data.append([int(row[0]), int(row[1])])

def extended_euclid_algorithm(a, b):
    assert a >= b and b >= 0
    if b == 0:
        return a, 1, 0
    d, x_curr, y_curr = extended_euclid_algorithm(b, a % b)
    x, y = y_curr, x_curr - y_curr * (a // b)
    return d, x, y

for index, d in enumerate(data):
    swap = d[0] < d[1]
    d[0], d[1] = max(d[0], d[1]), min(d[0], d[1])
    c, x, y= extended_euclid_algorithm(d[0], d[1])
    assert(x * d[0] + y * d[1] == c)
    if swap:
        x, y = y, x
    print(f"x={x},y={y},c={c}")