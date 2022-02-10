#include<iostream>
#include<vector>
#include<string>

const long long base = 1000;
const int base_len = 3;

class BigInteger {
private:
    bool sign = true;
    std::vector<long long> digits; // по основанию 1е9
public:
    BigInteger(long long n = 0); //конструктор от инт(неявный)
    explicit BigInteger(std::string s);//конструктор от строки
    void swap(BigInteger &n);
    bool Sign() const; //возвр знак
    size_t Size() const;
    void RemoveZeroes();
    long long &operator[](long long index);
    const long long &operator[](long long index) const;
    BigInteger operator-() const;
    BigInteger &operator+=(BigInteger n);
    BigInteger &operator-=(const BigInteger &n);
    BigInteger &operator*=(const BigInteger &n);
    BigInteger &operator/=(const BigInteger &n);
    BigInteger &operator%=(const BigInteger &n);
    BigInteger &operator++();//префиксный инкремент
    BigInteger operator++(int);//постфиксный
    BigInteger &operator--();//префиксный декремент
    BigInteger operator--(int);//постфиксный
    explicit operator bool() const; //неявное приведение к bool
    std::string toString() const;
    static BigInteger GCD(BigInteger n, BigInteger m);
};

BigInteger::BigInteger(long long n) { //check
    if(n == 0)
        digits.push_back(0);
    sign = (n >= 0);
    if(!sign) n *= -1;
    while (n > 0) {
        digits.push_back(n % base);
        n /= base;
    }
}

void BigInteger::swap(BigInteger &n) { //copy and swap
    std::swap(sign, n.sign);
    digits.swap( n.digits);
}

BigInteger::BigInteger(std::string s) {
   if(s.empty()) {
       *this = BigInteger();
       return;
   }
   sign = (s[0] != '-');
   if(s[0] == '-' || s[0] == '+') s.erase(0,1);
   while(s.size() > base_len) {
       digits.push_back(stoll(s.substr(s.size() - base_len)));
       s.erase(s.size() - base_len,s.size() - 1);
   }
   digits.push_back(stoll(s));
   RemoveZeroes();
}

bool BigInteger::Sign() const { //check
    return sign;
}

size_t BigInteger::Size() const { //check
    return digits.size();
}

void BigInteger::RemoveZeroes() {
    while (digits.size() > 1 && digits.back() == 0) digits.pop_back();
    if (digits.size() == 1 && digits[0] == 0) sign = true;
}

long long &BigInteger::operator[](long long index) { //вадратные скобки для неконст
    return digits[index];
}

const long long &BigInteger::operator[](long long index) const {//кв скобки для константных
    return digits[index];
}

bool operator<(const BigInteger &n, const BigInteger &m) {
   if(n.Sign() != m.Sign()) return(!n.Sign());
   if(n.Size() != m.Size()) return(n.Sign()) ^ (n.Size() > m.Size());
   for(long long i = static_cast<long long>(n.Size()) - 1; i >= 0; --i) {
       if(n[i] != m[i]) return (n.Sign() ^ (n[i] > m[i]));
   }
   return false;
}

bool operator>(const BigInteger &n, const BigInteger &m) {
    return m < n;
}

bool operator==(const BigInteger &n, const BigInteger &m) {
    return !(n < m || m < n);
}

bool operator!=(const BigInteger &n, const BigInteger &m) {
    return (n < m || m < n);
}

bool operator<=(const BigInteger &n, const BigInteger &m) {
    return (n < m || n == m);
}

bool operator>=(const BigInteger &n, const BigInteger &m) {
    return (m < n || n == m);
}

BigInteger BigInteger::operator-() const{
    BigInteger copy = *this;
    if(copy != 0)
        copy.sign = !copy.sign;
    return copy;
}

BigInteger &BigInteger::operator+=(BigInteger n) {
    int sgn = 1;
    if(sign != n.sign) {
        if(!sign) { //-*this + bi
            if(n > -*this) { //bi - *this > 0
                sign = true;
                swap(n);
            }
        }
        else { // *this + (-bi)
            if((*this) < -n) {
                swap(n);
            }
        }
        sgn = -1;
    }
    long long carry = 0;
    size_t i = 0;
    for(i = 0; i < n.Size() && i < Size(); ++i) {
        long long add = digits[i] + n.digits[i] * sgn + carry;
        digits[i] = add % base;
        carry = add / base;
        if(digits[i] < 0) {
            --carry;
            digits[i] += base;
        }
    }
    for(; i < Size() && carry != 0 ;++i) {
        long long add = digits[i] + carry;
        digits[i] = add % base;
        carry = add / base;
        if(add < 0) {
            digits[i] += base;
            --carry;
        }
    }
    for(; i < n.Size(); ++i) {
        long long add = n.digits[i] + carry;
        carry = add / base;
        digits.push_back(add % base);
    }
    if(carry != 0)
        digits.push_back(carry);
    RemoveZeroes();
    return *this;
}

BigInteger &BigInteger::operator-=(const BigInteger &n) {
    return (*this) += (-n);
}

BigInteger &BigInteger::operator*=(const BigInteger &n) {
    if(*this == 0 || n == 0)
        return *this = 0;
    size_t sz = digits.size() + n.Size();
    BigInteger res = 0;
    res.digits.clear();
    res.sign = (sign == n.sign);
    long long carry = 0;
    BigInteger small = *this, big = n;
    small.sign = true;
    big.sign = true;
    if(small > big) small.swap(big);
    for(size_t i = 0; i < sz + 1; ++i) {
        res.digits.push_back(carry);
        carry = 0;
        for (size_t j = 0; j <= i; ++j) {
            size_t a = 0, b = 0;
            if (j < big.digits.size())
                a = big.digits[j];
            if (i - j < small.digits.size())
                b = small.digits[i - j];
            res.digits[i] += a * b;
            carry += res.digits[i] / base;
            res.digits[i] %= base;
        }
    }
    *this = res;
    RemoveZeroes();
    return *this;
}

BigInteger &BigInteger::operator/=(const BigInteger &n) {
    BigInteger answer = 0;
    if (n == 0 || (*this) == 0)
        return *this = 0;
    BigInteger temp = n;
    temp.sign = true;
    BigInteger first = *this, sec = temp;
    first.sign = true;
    if(first < sec)
        return *this = 0;
    BigInteger prev_s = 1, prev_d = 1,  degree = 1;
    while (first >= sec) {
        degree = 1;
        prev_d = degree;
        prev_s = sec;
        while (first >= sec) {
            prev_s = sec;
            prev_d = degree;
            sec *= base;
            degree  *= base;
        }
        degree = prev_d;
        sec = prev_s;
        while (first >= sec) {
            first -= sec;
            answer += degree;
        }
        sec = temp;
        sec.RemoveZeroes();
    }
    answer.sign = (sign == n.sign);
    *this = answer;
    RemoveZeroes();
    return *this;
}

BigInteger &BigInteger::operator%=(const BigInteger &n) {
    BigInteger copy = *this;
    copy /= n;
    copy *= n;
    return *this -= copy;
}

BigInteger &BigInteger::operator++() {
    *this += 1;
    return *this;
}

BigInteger BigInteger::operator++(int) {
    BigInteger copy = *this;
    ++(*this);
    return copy;
}

BigInteger &BigInteger::operator--() {
    *this -= 1;
    return *this;
}

BigInteger BigInteger::operator--(int) {
    BigInteger copy = *this;
    --(*this);
    return copy;
}

BigInteger::operator bool() const {
    return(*this != 0);
}

std::string BigInteger::toString() const {
    std::string s;
    if (!sign) s = "-";
    else s = "";

    long long length = digits.size();
    s += std::to_string(digits.back());
    for (long long i = length - 2; i >= 0; --i) {
        std::string add = std::to_string(digits[i]);
        s += (std::string(base_len - add.size(), '0') + add);
    }
    return s;
}

BigInteger operator+(const BigInteger &n, const BigInteger &m) {
    BigInteger copy = n;
    copy += m;
    return copy;
}

BigInteger operator-(const BigInteger &n, const BigInteger &m) {
    BigInteger copy = n;
    copy -= m;
    return copy;
}

BigInteger operator*(const BigInteger &n, const BigInteger &m) {
    BigInteger copy = n;
    copy *= m;
    return copy;
}

BigInteger operator/(const BigInteger &n, const BigInteger &m) {
    BigInteger copy = n;
    copy /= m;
    return copy;
}

BigInteger operator%(const BigInteger &n, const BigInteger &m) {
    BigInteger copy = n;
    copy %= m;
    return copy;
}

std::ostream &operator<<(std::ostream &out, const BigInteger &n) { //вывод в поток
    return (out << n.toString());
}

std::istream &operator>>(std::istream &in, BigInteger &n) {
    std::string s;
    in >> s;
    n = BigInteger(s);
    return in;
}

BigInteger BigInteger::GCD(BigInteger n, BigInteger m) {
    n.sign = true;
    m.sign = true;
    if(n < m) n.swap(m);
    while(m) {
        n %= m;
        n.swap(m);
    }
    return n;
}

class Rational {
private:
    BigInteger numerator;
    BigInteger denominator;
    void Prime();
public:
    Rational(const BigInteger &n = 0, const BigInteger &d = 1);// КОНСТРУКТОР ОТ БИГ ИНТ
    Rational(long long n) : Rational(BigInteger(n)) {}; //конструктор от инт(неявный)

    std::string toString() const;

    explicit operator double() const;

    Rational operator-() const;

    Rational &operator+=(const Rational& r);

    Rational &operator-=(const Rational &r);

    Rational &operator*=(Rational r);

    Rational &operator/=(Rational r);

    std::string asDecimal(size_t precision = 0) const;

    friend bool operator<(const Rational &a, const Rational &b);

};

void Rational::Prime() {
    if(!denominator.Sign()) {
        denominator = -denominator;
        numerator = -numerator;
    }
    BigInteger x = BigInteger::GCD(numerator, denominator);
    numerator /= x;
    denominator /= x;
}

Rational::Rational(const BigInteger &n, const BigInteger &d) : numerator(n), denominator(d) {
    Prime();
}

std::string Rational::toString() const {
    std::string s = "";
    s += numerator.toString();
    if (denominator != BigInteger(1)) {
        s += "/";
        s += denominator.toString();
    }
    return s;
}

Rational Rational::operator-() const {
    return Rational(-numerator, denominator);
}

Rational &Rational::operator+=(const Rational& r) {
    numerator *= r.denominator;
    numerator += r.numerator * denominator;
    denominator *= r.denominator;
    Prime();
    return *this;
}

Rational operator+(const Rational &a, const Rational &b) {
    Rational copy = a;
    copy += b;
    return copy;
}

Rational &Rational::operator-=(const Rational &r) {
    *this += (-r);
    return *this;
}

Rational operator-(const Rational &a, const Rational &b) {
    Rational copy = a;
    copy -= b;
    return copy;
}

Rational &Rational::operator*=(Rational r) {
    denominator *= r.denominator;
    numerator *= r.numerator;
    Prime();
    return *this;
}

Rational operator*(const Rational &a, const Rational &b) {
    Rational copy = a;
    copy *= b;
    return copy;
}

Rational &Rational::operator/=(Rational r) {
    numerator *= r.denominator;
    denominator *= r.numerator;
    Prime();
    return *this;
}

Rational operator/(const Rational &a, const Rational &b) {
    Rational copy = a;
    copy /= b;
    return copy;
}

std::string Rational::asDecimal(size_t precision) const {
    BigInteger copy_num = numerator;
    BigInteger before_dot = copy_num / denominator;
    std::string ans = ((before_dot == 0 && !copy_num.Sign()) ? "-" : "") + before_dot.toString();
    if(!copy_num.Sign()) copy_num *= -1;
    copy_num %= denominator;
    if(precision > 0) {
        ans += ".";
        for(size_t i = 0; i < precision; ++i) {
            copy_num *= 10;
            ans += (copy_num / denominator).toString();
            copy_num %= denominator;
        }
    }
    return ans;
}

bool operator<(const Rational &a, const Rational &b) {
    return(a.numerator * b. denominator < a.denominator * b.numerator);
}

bool operator>(const Rational &a, const Rational &b) {
    return b < a;
}

bool operator==(const Rational &a, const Rational &b) {
    return !(b < a || a < b);
}

bool operator!=(const Rational &a, const Rational &b) {
    return (b < a || a < b);
}

bool operator<=(const Rational &a, const Rational &b) {
    return (a < b || a == b);
}

bool operator>=(const Rational &a, const Rational &b) {
    return (b < a || a == b);
}

Rational::operator double() const {
    double ans = std::stod(this->asDecimal(32));
    return ans;
}