import csv

def check_bound(x):
    return x >= 0 and x < 10**9

def extended_euclid_algorithm(a, b):
    assert a >= b and b >= 0
    if b == 0:
        return a, 1, 0
    d, x_curr, y_curr = extended_euclid_algorithm(b, a % b)
    x, y = y_curr, x_curr - y_curr * (a // b)
    return d, x, y

if __name__=='__main__':
    data = []

    with open('testinput-crt.txt', 'r') as f:
        reader = csv.reader(f, delimiter=',')
        for row in reader:
            data.append([int(row[0]), int(row[1]), int(row[2]), int(row[3])])

    for index, d in enumerate(data):
        a, m, b, n = d[0], d[1], d[2], d[3]
        assert check_bound(n) and check_bound(m) and check_bound(abs(b)) \
            and check_bound(abs(a)) and m != 0 and n != 0, "Bounds for m, n, a or b are not satisfied"

        '''
        To find x: x - a = m*y1 and x - b = n*y2 ==> my1 - ny2 = b - a
        By CRT: For gcd(m, n) = 1, there exists a unique solution % (m * n)
        '''
        if m < n:
            m, n = n, m
            a, b = b, a
        gcd_m_n, y1, y2_neg = extended_euclid_algorithm(m, n)
        if gcd_m_n != 1:
            print("-1")
        else:
            '''
            y1, y2_neg are solutions to my1 - ny2 = 1
            '''
            y1 *= (b - a)
            x = a + m * y1
            print(x % (m * n))


