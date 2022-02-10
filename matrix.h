#include<iostream>
#include<vector>

#include <iostream>
#include <vector>

#include<string>

using ll = long long;

const ll base = 1'000'000'000;
const int base_len = 9;
class BigInteger {
private:
    bool sign = true;
    std::vector<ll> digits; // по основанию 1е9
public:
    BigInteger(ll n = 0); //конструктор от инт(неявный)
    explicit BigInteger(std::string s);//конструктор от строки
    void swap(BigInteger &n);
    bool Sign() const; //возвр знак
    size_t Size() const;
    void RemoveZeroes();
    void Shift();
    ll &operator[](ll index);
    const ll &operator[](ll index) const;
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
    BigInteger Abs() const;
    void DivideByTwo();
    void MultiplyByTwo();
};

BigInteger::BigInteger(ll n) { //check
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

ll &BigInteger::operator[](ll index) { //вадратные скобки для неконст
    return digits[index];
}

const ll &BigInteger::operator[](ll index) const {//кв скобки для константных
    return digits[index];
}

bool operator<(const BigInteger &n, const BigInteger &m) {
    if(n.Sign() != m.Sign()) return(!n.Sign());
    if(n.Size() != m.Size()) return(n.Sign()) ^ (n.Size() > m.Size());
    for(ll i = static_cast<ll>(n.Size()) - 1; i >= 0; --i) {
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
    ll carry = 0;
    size_t i = 0;
    for(i = 0; i < n.Size() && i < Size(); ++i) {
        ll add = digits[i] + n.digits[i] * sgn + carry;
        digits[i] = add % base;
        carry = add / base;
        if(digits[i] < 0) {
            --carry;
            digits[i] += base;
        }
    }
    for(; i < Size() && carry != 0 ;++i) {
        ll add = digits[i] + carry;
        digits[i] = add % base;
        carry = add / base;
        if(add < 0) {
            digits[i] += base;
            --carry;
        }
    }
    for(; i < n.Size(); ++i) {
        ll add = n.digits[i] + carry;
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
    if(*this == 0 || n == 0) return *this = 0;
    if (n == base) {
        digits.insert(digits.begin(), 0);
        return *this;
    }

    if (n == BigInteger(-1)) {
        sign = !sign;
        return *this;
    }

    std::vector<ll> result(digits.size() + n.digits.size() + 1, 0);
    bool result_sign = (sign == n.sign);
    ll carry;
    for (size_t i = 0; i < digits.size(); i++) {
        carry = 0;
        size_t j = 0;
        while(j < n.digits.size()) {
            result[i + j] += carry + digits[i] * n.digits[j];
            carry = result[i + j] / base;
            result[i + j] %= base;
            ++j;
        }
        if (carry != 0)
            result[n.digits.size() + i] += carry;
    }
    digits = std::move(result);
    sign = result_sign;
    RemoveZeroes();
    return *this;
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

BigInteger &BigInteger::operator/=(const BigInteger &n) {//Бинпоиском
    if (Abs() < n.Abs())
        return (*this) = 0;
    BigInteger bi = n;
    bi.sign = true;
    BigInteger result, current;
    result.digits.resize(digits.size());
    for (ssize_t i = digits.size() - 1; i >= 0; --i) {
        if (current.digits.empty()) current.digits.push_back(0);
        if (!((current.digits.size() == 1) && (current.digits[0] == 0))) {
            current.digits.push_back(current.digits[current.digits.size() - 1]);
            for (ssize_t q = current.digits.size() - 2; q > 0; --q) current.digits[q] = current.digits[q - 1];
        }
        current.digits[0] = digits[i];
        current.RemoveZeroes();
        ll bi_sz = bi.digits.size(), cur_sz = current.digits.size(), x = 0;
        ll first_digit, second_digit;
        if(bi_sz < cur_sz) first_digit = current.digits[bi_sz];
        else first_digit = 0;
        if(bi_sz - 1< cur_sz) second_digit = current.digits[bi_sz - 1];
        else second_digit = 0;
        ll left = (base * first_digit+ second_digit) / (bi.digits.back() + 1);
        ll right = (base * first_digit + second_digit + 1) / bi.digits.back();
        while (left <= right) {
            ll middle = (left + right) / 2;
            BigInteger temp = bi * middle;
            if (temp <= current) {
                x = middle;
                left = middle  + 1;
            }
            else right = middle - 1;
        }
        result.digits[i] = x;
        current = current - bi * x;
    }
    result.RemoveZeroes();
    result.sign = (sign == n.sign);
    *this = result;
    return *this;
}

BigInteger& BigInteger::operator%=(const BigInteger& n) {
    if (Abs() < n.Abs())
        return *this;
    BigInteger bi = n;
    bi.sign = true;
    BigInteger result, current;
    result.digits.resize(digits.size());
    for (ssize_t i = digits.size() - 1; i >= 0; --i) {
        if (current.digits.empty()) current.digits.push_back(0);
        if (!((current.digits.size() == 1) && (current.digits[0] == 0))) {
            current.digits.push_back(current.digits[current.digits.size() - 1]);
            for (ssize_t q = current.digits.size() - 2; q > 0; --q) current.digits[q] = current.digits[q - 1];
        }
        current.digits[0] = digits[i];
        current.RemoveZeroes();
        ll bi_sz = bi.digits.size(), cur_sz = current.digits.size(), x = 0;
        ll first_digit, second_digit;
        if(bi_sz < cur_sz) first_digit = current.digits[bi_sz];
        else first_digit = 0;
        if(bi_sz - 1< cur_sz) second_digit = current.digits[bi_sz - 1];
        else second_digit = 0;
        ll left = (base * first_digit + second_digit) / (bi.digits.back() + 1);
        ll right = (base * first_digit + second_digit + 1) / bi.digits.back();
        while (left <= right) {
            ll middle = (left + right) / 2;
            BigInteger temp = bi * middle;
            if (temp <= current) {
                x = middle;
                left = middle  + 1;
            }
            else right = middle - 1;
        }
        result.digits[i] = x;
        current = current - bi * x;
    }
    current.RemoveZeroes();
    current.sign = sign;
    *this = current;
    return *this;
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

    ll length = digits.size();
    s += std::to_string(digits.back());
    for (ll i = length - 2; i >= 0; --i) {
        std::string add = std::to_string(digits[i]);
        s += (std::string(base_len - add.size(), '0') + add);
    }
    return s;
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
    if (n < m) std::swap(n, m);
    while (m > 0) {
        n %= m;
        std::swap(n, m);
    }
    return n;
}

void BigInteger::DivideByTwo() {
    if ((digits[0] % 2) == 1) return;
    for (size_t i = 0; i < digits.size(); ++i) {
        if ((digits[i] % 2) == 1) {
            digits[i-1] += base / 2;
        }
        digits[i] /= 2;
    }
    RemoveZeroes();
}

void BigInteger::MultiplyByTwo() {
    ll carry = 0;
    for (auto &digit : digits) {
        digit = digit * 2 + carry;
        if (digit >= base) {
            carry = 1;
            digit %= base;
        }
        else carry = 0;
    }
    if (carry == 1) digits.push_back(1);
}

BigInteger BigInteger::Abs() const {
    BigInteger copy = *this;
    copy.sign = true;
    return copy;
}



class Rational {
private:
    BigInteger numerator;
    BigInteger denominator;
    void Prime();
public:
    Rational(const BigInteger &n = 0, const BigInteger &d = 1);// КОНСТРУКТОР ОТ БИГ ИНТ
    Rational(ll n) : Rational(BigInteger(n)) {}; //конструктор от инт(неявный)

    Rational(std::string s);

    std::string toString() const;

    explicit operator double() const;

    Rational operator-() const;

    Rational &operator+=(const Rational& r);

    Rational &operator-=(const Rational &r);

    Rational &operator*=(const Rational &r);

    Rational &operator/=(const Rational &r);

    std::string asDecimal(size_t precision = 0) const;

    Rational abs();

    friend bool operator<(const Rational &a, const Rational &b);

};

Rational Rational::abs() {
    Rational copy = *this;
    if(!copy.denominator.Sign())
        copy.denominator = -copy.denominator;
    if(!copy.numerator.Sign())
        copy.numerator = -copy.numerator;
    return copy;
}

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

Rational::Rational(std::string s) {
    bool sgn = true;
    if (s[0] == '-') {
        sgn = false;
        s = s.substr(1, s.size() - 1);
    }
    ll slash = -1;
    for (size_t i = 0; i < s.size(); ++i)
        if (s[i] == '/') {
            slash = i;
            break;
        }
    if (slash == -1) {
        numerator = BigInteger(s);
        denominator = 1;
    } else {
        numerator = BigInteger(s.substr(0, slash));
        size_t length = s.size() - slash;
        denominator = BigInteger(s.substr(slash + 1, length));
    }
    if(!sgn)
        numerator = -numerator;
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
    numerator *= r.denominator;
    numerator -= r.numerator * denominator;
    denominator *= r.denominator;
    Prime();
    return *this;
}

Rational operator-(const Rational &a, const Rational &b) {
    Rational copy = a;
    copy -= b;
    return copy;
}

Rational &Rational::operator*=(const Rational &r) {
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

Rational &Rational::operator/=(const Rational &r) {
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

std::istream &operator>>(std::istream &in, Rational &r) {
    std::string s;
    in >> s;
    r = Rational(s);
    return in;
}


using std::vector;

template <int N>
class Finite {
private:
    ll elem;
    void MakeElemValid();
public:
    Finite() = default;
    Finite(ll _elem);
    Finite& operator=(const Finite &other);
    Finite &operator++();
    Finite operator-();
    bool operator==(const Finite& x) const;
    bool operator!=(const Finite& x) const;
    Finite &operator+=(const Finite& x);
    Finite &operator-=(const Finite& x);
    Finite &operator*=(const Finite& x);
    Finite &operator/=(const Finite& x);
    Finite pow(ll p) const;

    std::string toString() const;
};

template<int N>
Finite<N>& Finite<N>::operator=(const Finite<N> &other) {
    elem = other.elem;
    return *this;
}

template<int N>
void Finite<N>::MakeElemValid() {
    if(elem >= 0) elem %= N;
    else {
        elem = N - std::abs(elem) % N;
    }
}

template<int N>
Finite<N>::Finite(const ll _elem) {
    elem = _elem;
    MakeElemValid();
}

template<int N>
Finite<N>& Finite<N>::operator++() {
    ++elem;
    MakeElemValid();
    return *this;
}

template<int N>
bool Finite<N>::operator==(const Finite<N>& x) const{
    return elem == x.elem;
}

template<int N>
bool Finite<N>::operator!=(const Finite<N>& x) const{
    return !(elem == x.elem);
}

template<int N>
Finite<N>& Finite<N>::operator+=(const Finite<N>& x) {
    elem += x.elem;
    MakeElemValid();
    return *this;
}

template<int N>
Finite<N> operator+(const Finite<N> &first, const Finite<N> &second) {
    Finite<N> copy = first;
    copy += second;
    return copy;
}

template<int N>
Finite<N>& Finite<N>::operator-=(const Finite<N>& x) {
    elem -= x.elem;
    MakeElemValid();
    return *this;
}

template<int N>
Finite<N> operator-(const Finite<N> &first, const Finite<N> &second) {
    Finite<N> copy = first;
    copy -= second;
    return copy;
}

template<int N>
Finite<N> Finite<N>::operator-() {
    Finite<N> copy = Finite(0) - (*this);
    return copy;
}

template<int N>
Finite<N>& Finite<N>::operator*=(const Finite<N>& x) {
    elem *= x.elem;
    MakeElemValid();
    return *this;
}

template<int N>
Finite<N> operator*(Finite<N> first, const Finite<N> &second) {
    return first *= second;
}

template<int N>
Finite<N> Finite<N>::pow(ll p) const {
    if(p == 0)
        return Finite<N>(1);
    else if(p % 2 != 0)
        return (*this) * this->pow(p - 1);
    else {
        Finite<N> help = this->pow(p / 2);
        return help * help;
    }
}


template<int N> //TODO
Finite<N>& Finite<N>::operator/=(const Finite<N>& x) {
    //assert(isPrime(N));
    Finite<N> copy = x.pow(N - 2);
    elem *= copy.elem;
    return *this;
}

template<int N>
std::string Finite<N>::toString() const {
    return std::to_string(elem);
}

template<int N>
Finite<N> operator/(const Finite<N> &first, const Finite<N> &second) {
    Finite<N> copy = first;
    copy /= second;
    return copy;
}









// MATRIX

template <size_t M, size_t N, typename Field = Rational>
class Matrix {
protected:
    std::vector<std::vector<Field>> elements;
public:
    Matrix() {
        elements = std::vector<std::vector<Field>>(M, std::vector<Field>(N, static_cast<Field>(0)));
        if(M == N) {
            for(size_t i = 0; i < N; ++i)
                elements[i][i] = Field(1);
        }
    }
    Matrix(const std::vector<std::vector<Field>>& data);
    Matrix(const std::vector<std::vector<int>>& data);
    Matrix& operator=(const Matrix &other);
    bool operator==(const Matrix& another) const;
    bool operator!=(const Matrix& another) const;
    Matrix &operator+=(const Matrix& another);
    Matrix &operator-=(const Matrix& another);
    Matrix &operator*=(Field k);
    Matrix<M, M, Field> &operator*=(const Matrix<M, M, Field>& another);
    std::vector<Field>& operator[](size_t ind);
    const std::vector<Field>& operator[](size_t ind) const;
    std::vector<Field>& getRow(size_t ind);
    std::vector<Field> getColumn(size_t ind);
    void SwapRows(size_t from, size_t to);
    Matrix<N,M,Field> transposed() const;
    std::pair<Matrix<M, 2 * N, Field>, int>  Gauss() const;
    Field trace() const;
    Field det() const;
    size_t rank() const;
    Matrix& invert();
    Matrix inverted();
};

template<size_t M, size_t N, typename Field>
Matrix<M, N, Field>::Matrix(const std::vector<std::vector<Field>>& data) {
    elements.resize(M, std::vector<Field> (N, Field(0)));
    for (size_t row = 0; row < M; ++row)
        for (size_t column = 0; column < N; ++column)
            elements[row][column] = data[row][column];
}

template<size_t M, size_t N, typename Field >
Matrix<M, N, Field>::Matrix(const std::vector<std::vector<int>>& data) {
    elements.resize(M, std::vector<Field> (N, 0));
    for (size_t row = 0; row < M; ++row)
        for (size_t column = 0; column < N; ++column)
            elements[row][column] = data[row][column];
}

template<size_t M, size_t N, typename Field>
Matrix<M, N, Field>& Matrix<M, N, Field>::operator=(const Matrix &other) {
    for (size_t i = 0; i < M; ++i)
        for (size_t j = 0; j < N; ++j)
            elements[i][j] = other[i][j];
    return *this;
}


template<size_t M, size_t N, typename Field>
bool Matrix<M, N, Field>::operator==(const Matrix<M, N, Field> &another) const {
    Matrix<M, N, Field> copy = *this;
    for(size_t row = 0; row < M; ++row)
        for(size_t column = 0; column < N; ++column)
            if(elements[row][column] != another.elements[row][column])
                return false;
    return true;
}

template<size_t M, size_t N, typename Field>
bool Matrix<M, N, Field>::operator!=(const Matrix<M, N, Field> &another) const {
    return !(*this == another);
}

template<size_t M, size_t N, typename Field>
Matrix<M, N, Field>& Matrix<M, N, Field>::operator+=(const Matrix<M, N, Field> &another) {
    for(size_t row = 0; row < M; ++row)
        for(size_t column = 0; column < N; ++column)
            elements[row][column] += another.elements[row][column];
    return *this;
}

template<size_t M, size_t N, typename Field>
Matrix<M, N, Field>& operator+(const Matrix<M, N, Field>& first, const Matrix<M, N, Field> &another) {
    Matrix<M, N, Field> copy = first;
    return copy += another;
}

template<size_t M, size_t N, typename Field>
Matrix<M, N, Field>& Matrix<M, N, Field>::operator-=(const Matrix<M, N, Field> &another) {
    for(size_t row = 0; row < M; ++row)
        for(size_t column = 0; column < N; ++column)
            elements[row][column] -= another.elements[row][column];
    return *this;
}

template<size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator-(const Matrix<M, N, Field>& first, const Matrix<M, N, Field> &another) {
    Matrix<M, N, Field> copy = first;
    return copy -= another;
}

template<size_t M, size_t N, typename Field>
Matrix<M, N, Field>& Matrix<M, N, Field>::operator*=(const Field k) {
    for(size_t row = 0; row < M; ++row)
        for(size_t column = 0; column < N; ++column)
            elements[row][column] *= k;
    return *this;
}


template<size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator*(const Matrix<M, N, Field>& first, const Field k) {
    Matrix<M, N, Field> copy = first;
    return copy *= k;
}

template<size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator*(const Field k, const Matrix<M, N, Field>& first) {
    Matrix<M, N, Field> copy = first;
    return  copy *= k;
}

template<size_t M, size_t N, typename Field>
Matrix<M, M, Field>& Matrix<M, N, Field>::operator *= (const Matrix<M, M, Field>& another) {
    static_assert(M == N);
    std::vector< std::vector<Field> > help = std::vector< std::vector<Field> >(M, std::vector<Field>(M, 0));
    size_t ind = 0;
    while (ind < M) {
        for (size_t j = 0; j < M; ++j) {
            Field el = 0;
            for (size_t k = 0; k < M; ++k)
                el += elements[ind][k] * another[k][j];
            help[ind][j] = el;
        }
        ind++;
    }
    elements = help;
    return *this;
}


template<size_t M, size_t N, typename Field>
std::vector<Field>& Matrix<M, N, Field>::operator[](size_t ind) {
    return elements[ind];
}

template<size_t M, size_t N, typename Field>
const std::vector<Field>& Matrix<M, N, Field>::operator[](size_t ind) const{
    return elements[ind];
}

template<size_t M, size_t N, typename Field>
void Matrix<M, N, Field>::SwapRows(size_t from, size_t to) {
    std::swap(elements[from], elements[to]);
}

template<size_t M, size_t N, size_t K, typename Field = Rational>
Matrix<M,K,Field> operator*(const Matrix<M, N, Field> a, const Matrix<N, K, Field> b) {
    Matrix<M, K, Field> copy(std::vector<std::vector<Field>>(M, std::vector<Field> (K, 0)));
    for(size_t row = 0; row < M; ++row)
        for(size_t column = 0; column < K; ++column)
            for(size_t i = 0; i < N; ++i)
                copy[row][column] += a[row][i] * b[i][column];
    return copy;
}


template<size_t M, size_t N, typename Field>
std::vector<Field>& Matrix<M, N, Field>::getRow(size_t ind) {
    return elements[ind];
}

template<size_t M, size_t N, typename Field>
std::vector<Field> Matrix<M, N, Field>::getColumn(size_t ind) {
    std::vector<Field> answer;
    for(size_t i = 0; i < M; ++i)
        answer.push_back(elements[i][ind]);
    return answer;
}

template<size_t M, size_t N, typename Field>
Matrix<N, M, Field> Matrix<M, N, Field>::transposed() const{
    Matrix<N, M, Field> copy;
    for(size_t i = 0; i < M; ++i)
        for(size_t j = 0; j < N; ++j)
            copy[j][i] = elements[i][j];
    return copy;
}



template<size_t M, size_t N, typename Field>
std::pair<Matrix<M, 2 * N, Field>, int> Matrix<M, N, Field>::Gauss() const {
    int sign_of_det = 1;
    Matrix<M, 2 * N, Field> doubled_matrix;
    for(size_t row = 0; row < M; ++row) {
        for(size_t column = 0; column < N; ++column)
            doubled_matrix[row][column] = elements[row][column];
        doubled_matrix[row][row + N] = 1;
    }
    size_t ind = 0;
    for (size_t i = 0; (i < std::min(M, N)) && (ind < N); ++i) {
        for (size_t j = 0; j < M - i; ++j)
            if (doubled_matrix[i + j][ind] != 0) {
                if (j != 0) {
                    doubled_matrix.SwapRows(i, i + j);
                    sign_of_det = sign_of_det * (-1);
                }
                break;
            }
        if (doubled_matrix[i][ind] != 0) {
            for (size_t j = 1; j < M - i; ++j) {
                Field coefficient = (doubled_matrix[i + j][ind] / doubled_matrix[i][ind]);
                for (size_t k = 0; k < N - ind; ++k) {
                    doubled_matrix[i][ind + k] * coefficient;
                    doubled_matrix[i + j][ind + k] -= (doubled_matrix[i][ind + k] * coefficient);
                }
                for (size_t k = 0; k < N; ++k)
                    doubled_matrix[i + j][k + N] -= (doubled_matrix[i][k + N] * coefficient);
            }
            ++ind;
        }
    }
    return {doubled_matrix, sign_of_det};
}





template<size_t M, size_t N, typename Field>
Field Matrix<M, N, Field>::trace() const{
    static_assert(N == M);
    Field ans = 0;
    for(size_t i = 0; i < N; ++i)
        ans += elements[i][i];
    return ans;
}

template<size_t M, size_t N, typename Field>
Field Matrix<M, N, Field>::det() const{
    static_assert(N == M);
    std::pair<Matrix<N, 2 * N, Field>, int> ready = Gauss();
    Field det = 1;
    for (size_t i = 0; i < N; ++i)
        det *= ready.first[i][i];
    if (ready.second > 0)
        return det;
    return -det;
}

template<size_t M, size_t N, typename Field>
size_t Matrix<M, N, Field>::rank() const{
    std::pair<Matrix<M, 2 * N, Field>, int> ready = Gauss();
    size_t i = M - 1;
    bool flag = false;
    for (; i > 0; --i) {
        for (size_t j = 0; j < N; ++j)
            if (ready.first[i][j] != 0) {
                flag = true;
                break;
            }
        if (flag) break;
    }
    if (i != 0) return i + 1;
    for (size_t j = 0; j < N; ++j) if (ready.first[0][j] != 0) return 1;
    return 0;
}

template<size_t M, size_t N, typename Field>
Matrix<M, N, Field>& Matrix<M, N, Field>::invert() {
    static_assert(N == M);
    std::pair<Matrix<N, 2 * N, Field>, int> ready = Gauss();
    //assert(ready.first.MultiplyDia() != 0);

    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < i; ++j) {
            Field multiply = (ready.first[j][i] / ready.first[i][i]);
            for (size_t k = 0; k < N - i; ++k) {
                ready.first[j][i + k] -= ready.first[i][i + k] * multiply;
            }
            for (size_t k = 0; k < N; ++k) {
                ready.first[j][k + N] -= ready.first[i][k + N] * multiply;
            }
        }
        Field coefficient = Field(1) / ready.first[i][i];
        for (size_t j = i; j < 2 * N; ++j)
            ready.first[i][j] *= coefficient;
    }

    std::vector<std::vector<Field> > help = std::vector<std::vector<Field> >(N, std::vector<Field>(N, 0));
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < N; ++j)
            help[i][j] = ready.first[i][j + N];
    elements = help;
    return *this;
}

template<size_t M, size_t N, typename Field>
Matrix<M, N, Field> Matrix<M, N, Field>::inverted() {
    static_assert(N == M);
    Matrix<M, N, Field> copy = *this;
    return copy.invert();
}


template<size_t N, typename Field = Rational>
class SquareMatrix : public Matrix<N, N, Field> {
public:
    SquareMatrix(): Matrix<N, N, Field>() {
        Matrix<N, N , Field>::elements = std::vector<std::vector<Field>>(N, std::vector<Field>(N, 0));
        for (size_t i = 0; i < N; ++i)
            Matrix<N, N , Field>::elements[i][i] = Field(1);
    }

    explicit SquareMatrix(const std::vector<std::vector<Field>>& data): Matrix<N, N, Field>(data) {}

    explicit SquareMatrix(const std::vector<std::vector<int>>& data): Matrix<N, N, Field>(data) {}

};
