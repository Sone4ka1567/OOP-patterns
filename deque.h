#include <iostream>
#include <vector>
#include <cmath>

template <typename T>
class Deque {
private:
    std::vector<T*> arr;
    static const size_t inner_size = 64;
    std::pair<size_t, size_t> first_elem_ind = std::make_pair(0, 0); // первое значение - номер блока, второе - номер в блоке
    std::pair<size_t, size_t> last_elem_ind = std::make_pair(0, 0);
    size_t sz = 0;
public:

    Deque() {
        arr.resize(inner_size);
        first_elem_ind = std::make_pair(arr.size() / 2, inner_size / 2);
        last_elem_ind = std::make_pair(arr.size() / 2, inner_size / 2 - 1);
        for(size_t i = 0; i < arr.size(); ++i)
            arr[i] = reinterpret_cast<T*>(new int8_t[inner_size * sizeof(T)]);
    }

    Deque(const Deque& d): Deque() {
        first_elem_ind = d.first_elem_ind;
        last_elem_ind = d.last_elem_ind;
        sz = d.sz;
        arr.resize(d.arr.size());
        for(size_t i = 0; i < arr.size(); ++i) {
            arr[i] = reinterpret_cast<T*>(new int8_t[inner_size * sizeof(T)]);
            if(i >= first_elem_ind.first && i <= last_elem_ind.first) {
                for(size_t j = 0; j < inner_size; ++j) {
                    if(i == first_elem_ind.first && j >= first_elem_ind.second)  new(arr[i] + j) T(d.arr[i][j]);
                    else if (i == last_elem_ind.first && j <= last_elem_ind.second)  new(arr[i] + j) T(d.arr[i][j]);
                    else  new(arr[i] + j) T(d.arr[i][j]);
                }
            }
        }
    }

    explicit Deque(const int& n,  const T& value = T()) {
          size_t blocksCount = std::max(static_cast<size_t>(ceil(static_cast<double>(n) / static_cast<double>(inner_size))), static_cast<size_t>(1));
          arr.resize(blocksCount);
          first_elem_ind = std::make_pair(0, 0);
          last_elem_ind = std::make_pair(blocksCount - 1, (n - 1) % inner_size);
          for(size_t y = 0; y < blocksCount; ++y) {
              arr[y] = reinterpret_cast<T*>(new int8_t[inner_size * sizeof(T)]);
              for (size_t x = 0; x < inner_size; ++x) {
                  if(y == blocksCount - 1 && x > last_elem_ind.second) break;
                  new(arr[y] + x) T(value);
              }
          }
          sz = n;
    }

    Deque& operator=(const Deque& d) {
        Deque copy = Deque(d);
        swap(copy);
        return *this;
    }

    void swap(Deque& d) {
        std::swap(sz, d.sz);
        std::swap(arr, d.arr);
        std::swap(first_elem_ind, d.first_elem_ind);
        std::swap(last_elem_ind, d.last_elem_ind);
    }

    size_t size() const { return sz;}

    T& operator[](size_t i) {
       size_t pos_in_arr = first_elem_ind.first * inner_size + first_elem_ind.second + i;
       return arr[pos_in_arr / inner_size][pos_in_arr % inner_size];
    }

    const T& operator[](size_t i) const {
        size_t pos_in_arr = first_elem_ind.first * inner_size + first_elem_ind.second + i;
        return arr[pos_in_arr / inner_size][pos_in_arr % inner_size];
    }

    T& at(size_t i) {
        if(i >= sz) throw std::out_of_range("Error: out of range");
        size_t pos_in_arr = first_elem_ind.first * inner_size + first_elem_ind.second + i;
        return arr[pos_in_arr / inner_size][pos_in_arr % inner_size];
    }

    const T& at(size_t i) const{
        if(i >= sz) throw std::out_of_range("Error: out of range");
        size_t pos_in_arr = first_elem_ind.first * inner_size + first_elem_ind.second + i;
        return arr[pos_in_arr / inner_size][pos_in_arr % inner_size];
    }

    void push_back(const T& value) {
        if(last_elem_ind == std::make_pair(arr.size() - 1, inner_size - 1)) {
            size_t prev_arr_size = arr.size();
            arr.resize(2 * arr.size());
            for(size_t i = prev_arr_size; i < arr.size(); ++i) {
                arr[i] = reinterpret_cast<T*>(new int8_t[inner_size * sizeof(T)]);
            }
        }

        if(last_elem_ind.second == inner_size - 1) {
            last_elem_ind.first++;
            last_elem_ind.second = 0;
        }
        else last_elem_ind.second++;
        new(arr[last_elem_ind.first] + last_elem_ind.second) T(value);
        ++sz;
    }

    void pop_back() {
        (arr[last_elem_ind.first] + last_elem_ind.second)->~T();
        if(last_elem_ind.second == 0) {
            last_elem_ind.first--;
            last_elem_ind.second = inner_size - 1;
        }
        else last_elem_ind.second--;
        --sz;
    }

    void push_front(const T& value) {
        if(first_elem_ind == std::make_pair(size_t(0), size_t(0))) {
            size_t prev_size = arr.size();
            arr.resize(2 * prev_size);

            for(size_t i = 0; i < prev_size ; ++i) {
                arr[i + prev_size] = arr[i];
                arr[i] = reinterpret_cast<T*>(new int8_t[inner_size * sizeof(T)]);
            }
            first_elem_ind.first = prev_size;
            last_elem_ind.first += prev_size;
        }
        if(first_elem_ind.second == 0) {
            first_elem_ind.first--;
            first_elem_ind.second = inner_size - 1;
        }
        else first_elem_ind.second--;
        new(arr[first_elem_ind.first] + first_elem_ind.second) T(value);
        ++sz;
    }

    void pop_front() {
        (arr[first_elem_ind.first] + first_elem_ind.second)->~T();
        if(first_elem_ind.second == inner_size - 1) {
            first_elem_ind.first++;
            first_elem_ind.second = 0;
        }
        else first_elem_ind.second++;
        --sz;
    }

    //--------iterators------------//

    template <bool IsConst>
    struct common_iterator {
        typename std::vector<T*>::const_iterator ptr;
        size_t elem_ind;

        common_iterator(typename std::vector<T*>::const_iterator ptr, size_t elem_ind): ptr(ptr), elem_ind(elem_ind) {}

        common_iterator& operator++() {
            if(elem_ind == inner_size - 1) {
                ++ptr;
                elem_ind = 0;
            }
            else ++elem_ind;
            return *this;
        }

        common_iterator& operator--() {
            if(elem_ind == 0) {
                elem_ind = inner_size - 1;
                --ptr;
            }
            else --elem_ind;
            return *this;
        }

        int operator-(const common_iterator<IsConst>& second) const{
            return (ptr - second.ptr) * inner_size + static_cast<int>(elem_ind) - static_cast<int>(second.elem_ind);
        }

        common_iterator<IsConst> operator+(const int& n) const{
            int add_block = (elem_ind + n) / inner_size;
            int new_elem_index = (elem_ind + n) % inner_size;
            return common_iterator<IsConst>(ptr + add_block, new_elem_index);
        }

        common_iterator<IsConst> operator-(const int& n) const {
            int minus_block = (elem_ind - n) / inner_size;
            int new_elem_index = (elem_ind - n) % inner_size;
            if(new_elem_index < 0) {
                minus_block -= 1;
                new_elem_index += inner_size;
            }
            return common_iterator<IsConst>(ptr + minus_block, new_elem_index);
        }

        bool operator<(const common_iterator<IsConst>& second) const {
            return (*this - second) < 0;
        }

        bool operator>(const common_iterator<IsConst>& second) const {
            return (second - *this) < 0;
        }

        bool operator==(const common_iterator<IsConst>& second) const {
            return (!(*this < second) && !(*this > second));
        }

        bool operator!=(const common_iterator<IsConst>& second) const{
            return !(*this == second);
        }

        bool operator<=(const common_iterator<IsConst>& second) const{
            return !(*this > second);
        }

        bool operator>=(const common_iterator<IsConst>& second) const{
            return !(*this < second);
        }

        std::conditional_t<IsConst, const T&, T&> operator*() const {
            return *(*ptr + elem_ind);
        }

        std::conditional_t<IsConst, const T*, T*> operator->() const{
            return (*ptr + elem_ind);
        }

    };

    using iterator = common_iterator<false>;
    using const_iterator = const common_iterator<true>;

    iterator begin() {
        return iterator(arr.begin() + first_elem_ind.first, first_elem_ind.second);
    }

    iterator end() {
        auto copy = iterator(arr.begin() + last_elem_ind.first, last_elem_ind.second);
        ++copy;
        return copy;
    }

    const_iterator begin() const{
        return const_iterator(arr.begin() + first_elem_ind.first, first_elem_ind.second);
    }

    const_iterator end() const{
        auto copy = const_iterator(arr.begin() + last_elem_ind.first, last_elem_ind.second);
        return ++copy;
    }

    const_iterator cbegin() const {
        return const_iterator(arr.begin() + first_elem_ind.first, first_elem_ind.second);
    }

    const_iterator cend() const {
        auto copy = const_iterator(arr.begin() + last_elem_ind.first, last_elem_ind.second);
        return ++copy;
    }


    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    reverse_iterator rbegin() { return std::reverse_iterator(end());}

    reverse_iterator rend() { return std::reverse_iterator(begin());}

    const_reverse_iterator rbegin() const { return std::reverse_iterator(cend());}

    const_reverse_iterator rend() const { return std::reverse_iterator(cbegin());}

    const_reverse_iterator crbegin() const { return std::reverse_iterator(cend());}

    const_reverse_iterator crend() const { return std::reverse_iterator(cbegin());}

    void insert(iterator it, const T& value) {
        push_back(value);
        for(auto i = --end(); i != it; --i) {
            std::iter_swap(i, i - 1);
        }
    }

    void erase(iterator it) {
        for(auto i = it; i != --end(); ++i) {
            std::iter_swap(i, i + 1);
        }
        pop_back();
    }
};