import csv

def check_bound(x):
    return x > 0 and x < 10**9

def extended_euclid_algorithm(a, b):
    assert a >= b and b >= 0
    if b == 0:
        return a, 1, 0
    d, x_curr, y_curr = extended_euclid_algorithm(b, a % b)
    x, y = y_curr, x_curr - y_curr * (a // b)
    return d, x, y

def fast_modular_exponentiation(a, b, n):
    '''
    Returns: a^b modulo n
    '''
    if b == 0:
        return 1 % n
    
    half_power = fast_modular_exponentiation(a, b // 2, n)
    if(b % 2 == 0):
        return (half_power * half_power) % n
    else:
        return (((half_power * half_power) % n) * a) % n

if __name__=='__main__':
    data = []

    with open('testinput-Zn.txt', 'r') as f:
        reader = csv.reader(f, delimiter=',')
        for row in reader:
            data.append([int(row[0]), int(row[1]), int(row[2])])

    for index, d in enumerate(data):
        n, a, b = d[0], d[1], d[2]
        assert a < n and check_bound(n) and check_bound(a) and check_bound(b), "Bounds for n, a or b are not satisfied"

        a_power_b_mod_n = fast_modular_exponentiation(a, b, n)
        print(a_power_b_mod_n, end='')

        # Find solution to ax + ny = 1, x = a^(-1)
        gcd, y, x = extended_euclid_algorithm(n, a)
        if gcd == 1:
            print(f',true,{x % n}')
        else:
            print(',false')

