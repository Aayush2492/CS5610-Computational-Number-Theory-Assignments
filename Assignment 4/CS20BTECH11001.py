import csv

poly1 = []
poly2 = []

with open('input-polygcd1.csv', 'r') as f:
    reader = csv.reader(f, delimiter=',')
    p = [int(x) for x in next(reader)][0]
    poly1 = [int(x) for x in next(reader)][1:]
    poly2 = [int(x) for x in next(reader)][1:]

class Polynomial:
    def __init__(self, coeffs, p):
        self.p = p
        self.coeffs = [coeff % p for coeff in coeffs]
    
    def __repr__(self) -> str:
        return " + ".join([f"{coeff}*x^{i}" for i, coeff in enumerate(self.coeffs) if coeff != 0][::-1])
    
    def degree(self):
        # Returns index of highest degree non-zero coefficient otherwise 0
        degree = len(self.coeffs) - 1
        while degree >= 0 and self.coeffs[degree] == 0:
            degree -= 1
        return degree if degree >= 0 else -1
    
    def get_leading_coeff(self):
        index = self.degree()
        return self.coeffs[index] if index >= 0 else 0
    
    def add(self, other):
        # Z_p polynomial addition
        max_len = max(len(self.coeffs), len(other.coeffs))
        result = [0] * max_len

        poly1 = self.coeffs
        poly2 = other.coeffs

        for i in range(max_len):
            self_coeff = poly1[i] if i < len(poly1) else 0
            other_coeff = poly2[i] if i < len(poly2) else 0
            result[i] = (self_coeff + other_coeff) % self.p

        return Polynomial(result, self.p)
    
    def sub(self, other):
        # self - other
        # Z_p polynomial subtraction
        max_len = max(len(self.coeffs), len(other.coeffs))
        result = [0] * max_len

        poly1 = self.coeffs
        poly2 = other.coeffs

        for i in range(max_len):
            self_coeff = poly1[i] if i < len(poly1) else 0
            other_coeff = poly2[i] if i < len(poly2) else 0
            result[i] = (self_coeff - other_coeff) % self.p

        return Polynomial(result, self.p)
    
    def mul(self, other):
        # Z_p polynomial multiplication
        result = [0] * (len(self.coeffs) + len(other.coeffs) - 1)
        for i in range(len(self.coeffs)):
            for j in range(len(other.coeffs)):
                result[i + j] += self.coeffs[i] * other.coeffs[j]
        for i in range(len(result)):
            result[i] %= self.p
        return Polynomial(result, self.p)

    def divmod(self, other):
        # returns (q, r) such that self = q * other + r
        # Z_p polynomial division

        if other.degree() == -1:
            raise ValueError("Division by zero")
        
        if self.degree() < other.degree():
            return Polynomial([0], self.p), self
        
        # if other.degree() == 0:
        #     return Polynomial([pow(other.get_leading_coeff(), -1, self.p) * coeff for coeff in self.coeffs], self.p), Polynomial([0], self.p)
        
        q = Polynomial([0] * (self.degree() - other.degree() + 1), self.p)
        r = Polynomial(self.coeffs, self.p)

        while r.degree() >= other.degree():
            lead_coeff_r = r.get_leading_coeff()
            lead_coeff_other = other.get_leading_coeff()

            current_q = Polynomial([0] * (self.degree() - other.degree() + 1), self.p)
            current_q.coeffs[r.degree() - other.degree()] = (lead_coeff_r * pow(lead_coeff_other, -1, self.p)) % self.p
            r = r.sub(other.mul(current_q))

            q = q.add(current_q)

        return q, r

def extended_euclid_algorithm(a: Polynomial, b: Polynomial, p: int):
    # Assume all coefficients in Z_p
    assert a.degree() >= b.degree()
    if b.degree() == -1:
        return a, Polynomial([1], p), Polynomial([0], p)
    
    assert b.degree() >= 0, a.__repr__() + "||" + b.__repr__()
    
    q, r = a.divmod(b)
    # print("q", q)
    # print("r", r)
    
    d, x_curr, y_curr = extended_euclid_algorithm(b, r, p)
    x = y_curr
    y = x_curr.sub(q.mul(y_curr))
    return d, x, y

def post_process(gcd, u, v, p):
    # making gcd monic
    leading_coeff = gcd.get_leading_coeff()
    inverse_leading_coeff = pow(leading_coeff, -1, p)
    if leading_coeff != 1:
        gcd = Polynomial([inverse_leading_coeff * coeff for coeff in gcd.coeffs], p)
        u = Polynomial([inverse_leading_coeff * coeff for coeff in u.coeffs], p)
        v = Polynomial([inverse_leading_coeff * coeff for coeff in v.coeffs], p)
    return gcd, u, v

swap = False
if len(poly1) < len(poly2):
    poly1, poly2 = poly2, poly1
    swap = True

a = Polynomial(poly1[::-1], p)
b = Polynomial(poly2[::-1], p)

print("Polynomial f:", a)
print("Polynomial g:", b)

gcd, u, v = extended_euclid_algorithm(a, b, p)
gcd, u, v = post_process(gcd, u, v, p)
if swap:
    u, v = v, u
print("GCD:", gcd)
print("u:", u)
print("v:", v)