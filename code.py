class polynomial_arithmetic:
    # allows for the creation of a polynomial in GF(2^m) where m <= 163
    # the polynomial is represented as a list of coefficients
    def __init__(self, coefficients, p = 2):
        assert p == 2, "p must be 2" # we can uncomment this later if we want to implement other p
        self.m = len(coefficients) - 1 # in this case m can be 163
        self.p = p # p is always 2 for GF(2^m)
        self.original = coefficients
        self.coefficients = [x % p for x in coefficients][::-1]
        while self.coefficients and self.coefficients[-1] == 0: self.coefficients.pop() # remove leading 0s
        if not self.coefficients: self.coefficients = [0]
        self.coefficients.reverse()

    def __mul__(self, other):
        assert len(other.coefficients) > 0, "other must be a polynomial"
        new_other = other.coefficients[::-1]
        while len(new_other) < len(self.coefficients):
            new_other.append(0)
        new_other.reverse()
        # coefficients are either 0 or 1, so we can use XOR to make the multiplication faster
        # this is the same as multiplying polynomials in GF(2^m)
        final = polynomial_arithmetic([0])
        new_self = self.coefficients
        while new_self:
            if new_self[-1] == 1:
                final += polynomial_arithmetic(new_other)
            new_other.append(0)
            new_self.pop()
        return final


    def __add__(self, other):
        new_other = other.coefficients
        new_self = self.coefficients
        if len(new_other) > len(new_self):
            new_self, new_other = new_other, new_self
        new_other.reverse()
        while len(new_other) < len(new_self): new_other.append(0)
        new_other.reverse()
        return polynomial_arithmetic([x^y for x, y in zip(new_self, new_other)])
    
    def degree(self):
        return len(self.coefficients) - 1
    
    def __mod__(self, other):
        assert len(other.coefficients) > 0, "other must be a polynomial"
        ans = polynomial_arithmetic(self.coefficients)
        while ans.degree() >= other.degree():
            diff = ans.degree() - other.degree()
            ans += polynomial_arithmetic(other.coefficients + [0]*diff)
        return ans
    
    def __floordiv__(self, other):
        assert len(other.coefficients) > 0, "other must be a polynomial"
        ans = polynomial_arithmetic([0])
        while self.degree() >= other.degree():
            diff = self.degree() - other.degree()
            ans += polynomial_arithmetic([1] + [0]*diff)
            self += polynomial_arithmetic(other.coefficients + [0]*diff)
        return ans

    def __sub__(self, other):
        return self + other
    
    def val(self, bin = False):
        if bin: return "".join([str(x) for x in self.coefficients])
        ans = ""
        for i, x in enumerate(self.coefficients[::-1]):
            if x == 0: continue
            if i == 0: ans += str(x)
            elif i == 1: ans += f" + x" if ans else "x"
            else: ans += f" + x^{i}"
        return ans if ans else "0"
