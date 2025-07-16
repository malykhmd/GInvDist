#1661371678625669441005231671474323137508345315870356275200000000000000000000   
import numpy
import math
import time

def first_primes(n):
    primes = [2, 3, 5, 7]
    x = 0
    for i in range(10, n):
        j = 0
        fl = False
        while (not fl) and j < len(primes) and primes[j]**2 < i:
            if i % primes[j] == 0:
                fl = True
                break
            else:
                j += 1
        if not fl:
            primes.append(i)
            x = (x + 1) % 100000
            if x == 0:
                print(primes[-1]/n)
            
    return primes

def calc_primes():
    primes = first_primes(10000000)
    numpy.save('primes.npy', primes)
    return primes

def load_primes():
    primes = numpy.load('primes.npy')
    return primes

def primfacs(n):
    i = 0
    primfac = []
    while i < primes_count and primes[i]**2 <= n:
        while n % int(primes[i]) == 0:
            primfac.append(int(primes[i]))
            n = n // int(primes[i])
        i = i + 1
    if n > 1:
        primfac.append(n)
    return primfac

class Coef:
    def __init__(self):
        self.factors = []
        self.sign = 1
        self.reminder = None
        
    def __str__(self):
        if self.sign == 0:
            return '0'
        line = '-' if self.sign == -1 else ''
        for f in self.factors:
            line += str(primes[f[0]]) + ('^' + str(f[1]) if f[1] > 1 else '') + '*'
        line = line[:-1]
        line += '' if self.reminder == None else '*' + str(self.reminder)
        return line
    
    def __repr__(self):
        return ('-' if self.sign == -1 else '') + str(self.factors) + \
                ('' if self.reminder == None else '  r{' + str(self.reminder) + '}')
                
    def get_int_value(self):
        number = self.sign
        for f in self.factors:
            number *= pow(int(primes[f[0]]), f[1])
        if self.reminder:
            number *= self.reminder        
        return number
    
    def parse_int(int_value):
        coef = Coef()
        if int_value > 0:
            coef.sign = 1
        elif int_value < 0:
            coef.sign = -1
            int_value = -int_value
        else:
            coef.sign = 0
            return coef
        i = 0
        pow = 0
        while i < primes_count and int_value > 1:
            if int_value % int(primes[i]) == 0:
                pow += 1
                int_value //= int(primes[i])
            else: 
                if pow > 0:
                    coef.factors.append([i,pow])
                    pow = 0
                i += 1
        if pow > 0:
            coef.factors.append([i,pow])
            pow = 0
        if int_value > 1:
            coef.reminder = int_value
        return coef
    
    def gcd(self, other):
        assert isinstance(other, Coef)
        result = Coef()
        i = j = 0
        while i < len(self.factors) and j < len(other.factors):
            if self.factors[i][0] == other.factors[j][0]:
                result.factors.append([self.factors[i][0], min(self.factors[i][1], other.factors[j][1])])
                i += 1
                j += 1
            elif self.factors[i][0] < other.factors[j][0]:
                i += 1
            else:
                j += 1
        return result
    
    def rem(self, other):
        assert isinstance(other, Coef)
        rem = Coef()
        rem.sign = self.sign
        rem.reminder = self.reminder
        i = j = 0
        while i < len(self.factors) and j < len(other.factors):
            if self.factors[i][0] == other.factors[j][0]:
                if self.factors[i][1] > other.factors[j][1]:
                    rem.factors.append([self.factors[i][0], self.factors[i][1] - other.factors[j][1]])
                i += 1
                j += 1
            elif self.factors[i][0] < other.factors[j][0]:
                rem.factors.append([self.factors[i][0], self.factors[i][1]])
                i += 1
            else:
                j += 1
        while i < len(self.factors):
            rem.factors.append([self.factors[i][0], self.factors[i][1]])
            i += 1
        return rem
        
    def gcd_rems(self, other):
        assert isinstance(other, Coef)
        g = self.gcd(other)
        return g, self.rem(g), other.rem(g)
    
    def __add__(self, other):
        assert isinstance(other, Coef)
        if self.sign == 0:
            return other
        if other.sign == 0:
            return self
        
        result = Coef()
        rem1 = Coef()
        rem1.sign = self.sign
        rem1.reminder = self.reminder
        rem2 = Coef()
        rem2.sign = other.sign
        rem2.reminder = other.reminder
        
        i = j = 0
        while i < len(self.factors) and j < len(other.factors):
            if self.factors[i][0] == other.factors[j][0]:
                if self.factors[i][1] == other.factors[j][1]:
                    result.factors.append([self.factors[i][0], self.factors[i][1] + other.factors[j][1]])
                if self.factors[i][1] > other.factors[j][1]:
                    result.factors.append([self.factors[i][0], other.factors[j][1]])
                    rem1.factors.append([self.factors[i][0], self.factors[i][1] - other.factors[j][1]])
                if self.factors[i][1] < other.factors[j][1]:
                    result.factors.append([self.factors[i][0], self.factors[i][1]])
                    rem2.factors.append([self.factors[i][0], other.factors[j][1] - self.factors[i][1]])
                i += 1
                j += 1
            elif self.factors[i][0] < other.factors[j][0]:
                rem1.factors.append([self.factors[i][0], self.factors[i][1]])
                i += 1
            else:
                rem2.factors.append([other.factors[j][0], other.factors[j][1]])
                j += 1
        
        while i < len(self.factors):
            rem1.factors.append([self.factors[i][0], self.factors[i][1]])
            i += 1
        
        while j < len(other.factors):
            rem2.factors.append([other.factors[j][0], other.factors[j][1]])
            j += 1
        
        rem_sum = Coef.parse_int(rem1.get_int_value() + rem2.get_int_value())
        return result * rem_sum
    
    def __sub__(self, other):
        assert isinstance(other, Coef)
        if self.sign == 0:
            other.sign = -other.sign
            return other
        if other.sign == 0:
            return self
        result, rem1, rem2 = self.gcd_rems(other)
        rem_sum = Coef.parse_int(rem1.get_int_value() - rem2.get_int_value())
        return result * rem_sum
    
    def __mul__(self, other):
        assert isinstance(other, Coef)
        result = Coef()
        result.sign = self.sign * other.sign
        if result.sign == 0:
            return result
        
        i = j = 0
        while i < len(self.factors) and j < len(other.factors):
            if self.factors[i][0] == other.factors[j][0]:
                result.factors.append([self.factors[i][0], self.factors[i][1] + other.factors[j][1]])
                i += 1
                j += 1
            elif self.factors[i][0] < other.factors[j][0]:
                result.factors.append([self.factors[i][0], self.factors[i][1]])
                i += 1
            else:
                result.factors.append([other.factors[j][0], other.factors[j][1]])
                j += 1
        
        while i < len(self.factors):
            result.factors.append([self.factors[i][0], self.factors[i][1]])
            i += 1
        
        while j < len(other.factors):
            result.factors.append([other.factors[j][0], other.factors[j][1]])
            j += 1
        
        if self.reminder and other.reminder:
            result.reminder = self.reminder * other.reminder
        elif self.reminder:
            result.reminder = self.reminder
        else:
            result.reminder = other.reminder
            
        return result
        
    def __eq__(self, other):
        assert isinstance(other, Coef)
        if (self.sign != other.sign) or (self.reminder != other.reminder) or (len(self.factors) != len(other.factors)):
            return False
        for i in range(len(self.factors)):
            if self.factors[i][0] != other.factors[i][0] or self.factors[i][1] != self.factors[i][1]:
                return False
        return True

    def __ne__(self, other):
        return not self == other
    
    
def test1():
    #number = 1661371678625669441005231671474323137508345315870356275200000000000000000000    
    number = 1661371678625669441005231671474323137508345315870356275200000000000000000000

    now = time.time()
    print('Факторизация ')
    print(primfacs(number))
    print(time.time() - now)

    now = time.time()
    C = Coef.parse_int(number)
    print('Создание коэф. ')
    print(C)
    print(time.time() - now)
    
    now = time.time()
    print('Вычисление значения ')
    print(C.get_int_value())
    print(time.time() - now)

def test2():
    #number = 1661371678625669441005231671474323137508345315870356275200000000000000000000    
    num1 = 694410052316714
    num2 = 453158703562
    
    c1 = Coef.parse_int(num1)
    c2 = Coef.parse_int(num2)
    
    print('Сложение')
    
    now = time.time()
    num3 = num1 + num2
    print('Время обычного ', time.time() - now)
    
    now = time.time()
    c3 = c1 + c2
    print('Время факторизованного ', time.time() - now)
    
    print(c1, c1.get_int_value(), '==', num1)
    print(c2, c2.get_int_value(), '==', num2)
    print(c3, c3.get_int_value(), '==', num3)
    assert c3.get_int_value() == num3
    
    print('Вычитание')
    
    now = time.time()
    num3 = num1 - num2
    print('Время обычного ', time.time() - now)
    
    now = time.time()
    c3 = c1 - c2
    print('Время факторизованного ', time.time() - now)

    print(c1, c1.get_int_value(), '==', num1)
    print(c2, c2.get_int_value(), '==', num2)
    print(c3, c3.get_int_value(), '==', num3)
    assert c3.get_int_value() == num3
    
    print('Умножение')
    
    now = time.time()
    num3 = num1 * num2
    print('Время обычного ', time.time() - now)
    
    now = time.time()
    c3 = c1 * c2
    print('Время факторизованного ', time.time() - now)
    
    print(c1, c1.get_int_value(), '==', num1)
    print(c2, c2.get_int_value(), '==', num2)
    print(c3, c3.get_int_value(), '==', num3)
    assert c3.get_int_value() == num3
    print('Прошло')
    
    print('GCD')
    
    now = time.time()
    num3 =  math.gcd(num1, num2)
    print('Время обычного ', time.time() - now)
    
    now = time.time()
    c3 = c1.gcd(c2)
    print('Время факторизованного ', time.time() - now)
    
    assert c3.get_int_value() == num3
    print('Прошло')

primes = load_primes()
primes_count = len(primes)

if __name__ == '__main__':
    test2()
    
    