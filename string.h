#include <iostream>
#include <cstring>
class String {
private:
    size_t sz = 0;
    size_t capacity = 0;
    char* str = nullptr;
public:
    String(const char* s): sz(strlen(s)) ,capacity(2 * strlen(s)), str(new char[++capacity]) { // конструктор от С-Style
        memcpy(str, s, sz);
    }

    String(size_t n, char c = '\0'): sz(n),capacity(2 * n), str(new char[++capacity]) { //конструктор от числа и символа
        memset(str, c, sz);
    }

    String() = default; //конструктор по умолчанию

    String(const String& s) : String(s.sz, '\0') { //конструктор копирования делегирующий
        memcpy(str, s.str, sz);
    }

    String& operator=(const String& s) { //оператор присваивания
        String copy = s;
       swap(copy);
       return *this;
    }

    void swap(String& s) { //copy and swap
        std::swap(sz, s.sz);
        std::swap(str, s.str);
        std::swap(capacity, s.capacity);
    }

    char& operator[](size_t index) { //вадратные скобки для неконст
        return str[index];
    }

    const char& operator[](size_t index) const { //кв скобки для константных
        return str[index];
    }

    size_t length() const { //константный метод, работает для константных и нет
        return sz; //все поля считаются константными, менять ничего нельзя
    }

   void push_back(char c) { //do again
        if(sz >= capacity) {
            String copy = *this;
            swap(copy);
        }
        str[sz] = c;
        sz++;
    }

    void pop_back() { //pop_back()
        sz--;
        if (sz * 4 < capacity) {
                String copy = *this;
                swap(copy);
        }
    }

    String& operator+=(char c) {//оператор += от символа
        push_back(c);
        return *this;
    }

    String& operator+=(const String& s) { //оператор += от строки
        if(sz + s.sz >= capacity) {
           String copy = String(sz + s.sz, '\0');
           memcpy(copy.str, str, sz);
           memcpy(copy.str + sz, s.str, s.sz);
           swap(copy);
        }
        else {
            memcpy(str + sz, s.str, s.sz);
            sz += s.sz;
        }
        return *this;
    }

    size_t find(const String& substring) const {
        for(size_t i = 0; i < sz - substring.sz + 1; ++i) {
            if(str[i] == substring[0]) {
                bool flag = true;
                for(size_t j = 1; j < substring.sz; ++j) {
                    if(str[i+j] != substring[j]) {
                        flag = false;
                        break;
                    }
                }
                if(flag) {
                    return i;
                }
            }
        }
        return sz;
    }

    size_t rfind(const String& substring) const {
        for(size_t i = sz - 1; i > substring.sz - 1; --i) {
            if(str[i] == substring[substring.sz - 1]) {
                bool flag = true;
                for(size_t j = substring.sz - 2; j > 0; --j) {
                    if(str[i - substring.sz + j + 1] != substring[j]) {
                        flag = false;
                        break;
                    }
                }
                if(substring[0]!= str[i - substring.sz + 1]) {
                    flag = false;
                }
                if(flag) {
                    return i - substring.sz + 1;
                }
            }
        }
        return sz;
    }

    String substr(size_t start, size_t count) const {
        String ans = String(count, '\0');
        memcpy(ans.str, str + start, count);
        return ans;
    }

    char& front() { //метод front()
       return str[0];
   }

    const char& front() const{ //метод front() для констант
        return str[0];
    }

    const char& back() const{ // метод back() для констант
        return str[sz - 1];
    }

    char& back() { // метод back()
       return str[sz - 1];
   }


   bool empty() const {//empty()
       return sz == 0;
   }


   void clear() { //clear()
       String copy;
       swap(copy);
   }

    ~String() { //деструктор
        delete[] str;
    };

};

bool operator==(const String& s1, const String& s2) { //оператор ==
    if(s1.length() != s2.length())
        return false;
    for(size_t i = 0; i < s1.length(); ++i) {
        if(s1[i] != s2[i])
            return false;
    }
    return true;
}

String operator+(const String& s1, const String& s2) { //сложить строку со строкой
    //создает копию поэтому сначала +=
    String copy = s1; //оператор только правый обьект
    return copy += s2;
}

String operator+(const String& s, char c) { //сложить строку с символом
    String copy = s;
    return copy += c;
}

String operator+(char c, const String& s) { //сложить символ со строкой
    String copy = String(1, c);
    return copy + s;
}

std::ostream& operator<<(std::ostream& out, const String& s) { //вывод в поток
    for (size_t i = 0; i < s.length(); ++i) {
        out << s[i];
    }
    return out;
}

std::istream& operator>>(std::istream& in, String& s) {
    while (strchr(" \n\r\t", in.peek()) != nullptr) in.get();
    for(char c; strchr(" \n\r\t", in.peek()) == nullptr && in >> c;) {
        s.push_back(c);
    }
    return in;
}
