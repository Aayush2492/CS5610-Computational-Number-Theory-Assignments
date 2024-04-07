import csv

def legendre_symbol(a, p):
    '''
    Returns the Legendre symbol (a/p) = a^((p-1)/2) mod p
    '''
    return pow(a, (p - 1) // 2, p)


def find_k(b, p):
    '''
    returns the least non-negative integer k such that b^(2^k) = 1 (mod p)
    '''
    k = 0
    while True:
        e = 1 << k
        if pow(b, e, p) == 1:
            break
        k += 1
    return k


def tonelli_shanks(a, p):
    '''
    Returns the square root of a modulo p if it exists, else returns 0
    Constraints: p is prime and (a > 0 and a < p)
    '''
    if p == 2:
        return 1

    if legendre_symbol(a, p) != 1:
        return 0
    
    # Find p - 1 = 2^t * m where m is odd
    m, t = p - 1, 0
    while m % 2 == 0:
        m //= 2
        t += 1
    
    b = pow(a, m, p)
    k = find_k(b, p)
    
    if k == t:
        return 0
    
    x = pow(a, (m + 1) // 2, p)
    
    if k == 0:
        return x
    
    r = 2
    while legendre_symbol(r, p) != p - 1:
        r += 1
    s = pow(r, m, p)
    S = pow(s, pow(2, t - k), p)

    while k > 0:
        assert (x * x) % p == (a * b) % p, "Invariant failed at beginning of loop"
        b = (b * S) % p
        x = (x * pow(s, pow(2, t - k - 1), p)) % p
        k = find_k(b, p)
        S = pow(s, pow(2, t - k), p)
        assert (x * x) % p == (a * b) % p, "Invariant failed at end of loop"
    
    return x


if __name__ == "__main__":
    data = []
    with open('inputSquareRoots.csv', 'r') as f:
        reader = csv.reader(f, delimiter=',')
        for row in reader:
            data.append([int(float(row[0])), int(float(row[1]))])

    for d in data:
        a, p = d[0], d[1]
        assert a < p and a > 0, "Invalid input"
        ans = tonelli_shanks(a, p)
        ans = min(ans, p - ans)
        if ans != 0:
            assert (ans * ans) % p == a, "Output does not satisfy x^2 = a (mod p)"
        print(ans)