#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <memory>
#include <iterator>
#include <type_traits>

//-----FixedAllocator-----//

template<int chunk_size>
class FixedAllocator {
private:

    struct Node {
        char data[chunk_size];
        Node* next;
    };
    Node* free_node = nullptr;
    const int pool_size = 128;

    void create_new_pool() {
        Node* new_pool = new Node[pool_size];
        for (int ind = 1; ind < pool_size; ++ind)
            new_pool[ind - 1].next = new_pool + ind;
        new_pool[pool_size - 1].next = free_node;
        free_node = new_pool;
    }

public:
    static FixedAllocator& get_instance() {
        static FixedAllocator instance;
        return instance;
    }

    FixedAllocator() = default;

    void* allocate() {
        if (free_node == nullptr)
            create_new_pool();
        void* result = static_cast<void*>(free_node);
        free_node = free_node->next;
        return result;
    }

    void deallocate(void* ptr) {
        auto new_node = static_cast<Node*>(ptr);
        new_node->next = free_node;
        free_node = new_node;
    }
    ~FixedAllocator() = default;
};

//-----FastAllocator-----//

template<typename T>
class FastAllocator {
private:

    T* alloc_with_fixed_allocator() {return static_cast<T*>(FixedAllocator<sizeof(T)>::get_instance().allocate());}

    T* alloc_with_std_allocator(int n) {return reinterpret_cast<T*>(new char[n * sizeof(T)]);}

public:
    using value_type = T;

    template<typename U>
    using rebind = FastAllocator<U>;

    FastAllocator() = default;

    template<typename U>
    FastAllocator(const FastAllocator<U>&) {}

    T* allocate(int n) {
        if (n == 1) return alloc_with_fixed_allocator();
        return alloc_with_std_allocator(n);
    }

    void deallocate(T* ptr, int n) {
        if (n == 1) return FixedAllocator<sizeof(T)>::get_instance().deallocate(ptr);
        delete ptr;
    }

    ~FastAllocator() = default;
};

//-----List-----//

template<typename T, typename Allocator = std::allocator<T>>
class List {
private:
    struct Node {
        T value;
        Node* prev;
        Node* next;

        explicit Node(const T& value, Node* prev, Node* next) :value(value), prev(prev), next(next) {};
        Node(Node* prev = nullptr, Node* next = nullptr): value(T()), prev(prev), next(next) {}
    };

    using NodeAlloc = typename std::allocator_traits<Allocator>::template rebind_alloc<Node>;

    Node* head;
    int sz;
    NodeAlloc node_alloc;

    void link_nodes(Node* left, Node* right) {
        if (left != nullptr)
            left->next = right;
        if (right != nullptr)
            right->prev = left;
    }

    void destroy_deallocate_node(Node* node) {
        std::allocator_traits<NodeAlloc>::destroy(node_alloc, node);
        std::allocator_traits<NodeAlloc>::deallocate(node_alloc, node, 1);
    }

    Node* get_head() {
        Node* node_ = static_cast<Node*>(std::allocator_traits<NodeAlloc>::allocate(node_alloc, 1));
        node_->prev = node_;
        node_->next = node_;
        return node_;
    }

    //-----Iterators-----//

    template<bool Const>
    class common_iterator {
    private:
        Node* node;
    public:

        using iterator_category = std::bidirectional_iterator_tag;
        using value_type = T;
        using difference_type = int64_t;
        using pointer = typename std::conditional_t<Const, T const*, T*>;
        using reference = typename std::conditional_t<Const, T const&, T&>;

        explicit common_iterator(const Node* node) : node(const_cast<Node *>(node)) {}
        common_iterator() : node(nullptr) {};
        common_iterator(const common_iterator<false>& iter) : node(iter.node) {};

        bool operator==(const common_iterator& other) const {
            if (node == nullptr || other.node == nullptr)
                return false;
            return node == other.node;
        }
        bool operator!=(const common_iterator& other) const {return !(*this == other);}
        typename common_iterator::reference operator*() const {return node->value;}
        typename common_iterator::pointer operator->() const {return &node->value;}
        common_iterator& operator++() {
            node = node->next;
            return *this;
        }
        common_iterator operator++(int) {
            auto result = *this;
            ++(*this);
            return result;
        }

        common_iterator& operator--() {
            node = node->prev;
            return *this;
        }
        friend class List<T, Allocator>;
    };

public:
    using iterator = common_iterator<false>;
    using const_iterator = common_iterator<true>;

    iterator begin() { return iterator(head->next); }
    const_iterator begin() const { return const_iterator(head->next); }
    const_iterator cbegin() const { return const_iterator(head->next); }
    iterator end() { return iterator(head); }
    const_iterator end() const { return const_iterator(head); }
    const_iterator cend() const { return const_iterator(head); }

    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    reverse_iterator rbegin() { return std::reverse_iterator(end()); }
    reverse_iterator rend() { return std::reverse_iterator(begin()); }
    const_reverse_iterator rbegin() const{ return std::reverse_iterator(cend()); }
    const_reverse_iterator rend() const{ return std::reverse_iterator(cbegin()); }
    const_reverse_iterator crbegin() const{ return std::reverse_iterator(cend()); }
    const_reverse_iterator crend() const{ return std::reverse_iterator(cbegin()); }



    //-----List Methods-----//

    explicit List(const Allocator& alloc = Allocator()) : head(get_head()), sz(0), node_alloc(alloc) {}

    List(int count, const T& value, const Allocator& alloc = Allocator()) : head(get_head()),
                                                                            sz(0), node_alloc(alloc) {
        for (int index = 0; index < count; ++index)
            insert_before(begin(), value);
    }

    explicit List(int count, const Allocator& alloc = Allocator()) : head(get_head()),
                                                                     sz(0), node_alloc(alloc) {
        for(int i = 0; i < count; ++i) {
            auto prev = std::prev(begin()).node;
            auto new_node = std::allocator_traits<NodeAlloc>::allocate(node_alloc, 1);
            std::allocator_traits<NodeAlloc>::construct(node_alloc, new_node, prev, begin().node);
            link_nodes(new_node, begin().node);
            link_nodes(prev, new_node);
            ++sz;
        }
    }

    List(const List& other) : head(get_head()), sz(0), node_alloc(std::allocator_traits<Allocator>::select_on_container_copy_construction(other.node_alloc)) {
        for (const auto& elem : other)
            this->push_back(elem);
    }

    List& operator=(const List& other) {
        if (this == &other) return *this;
        auto current_node = head->next;
        while (current_node != head) {
            auto next = current_node->next;
            destroy_deallocate_node(current_node);
            current_node = next;
        }
        head->prev = head;
        head->next = head;
        sz = 0;
        if (std::allocator_traits<NodeAlloc>::propagate_on_container_copy_assignment::value) node_alloc = other.node_alloc;
        for (const auto& elem : other) push_back(elem);
        return *this;
    }

    NodeAlloc get_allocator() {return node_alloc;}

    void insert_before(const_iterator pos, const T& value) {
        auto prev = std::prev(pos).node;
        auto new_node = std::allocator_traits<NodeAlloc>::allocate(node_alloc, 1);
        std::allocator_traits<NodeAlloc>::construct(node_alloc, new_node, value, prev, pos.node);
        link_nodes(new_node, pos.node);
        link_nodes(prev, new_node);
        ++sz;
    }

    void erase(const_iterator pos) {
        if (pos != end()) {
            auto next = std::next(pos).node;
            link_nodes(std::prev(pos).node, next);
            destroy_deallocate_node(pos.node);
            --sz;
        }
    }

    void push_back(const T& value) {insert_before(end(), value);}
    void push_front(const T& value) { insert_before(begin(), value);}
    void pop_back() {erase(std::prev(end()));}
    void pop_front() { erase(begin());}
    void insert(const_iterator pos, const T& value) { insert_before(pos, value);}

    int size() const { return sz; }

    ~List() {
        auto current_node = head->next;
        while (current_node != head) {
            auto next = current_node->next;
            destroy_deallocate_node(current_node);
            current_node = next;
        }
        head->prev = head->next = head;
        sz = 0;
        std::allocator_traits<NodeAlloc>::deallocate(node_alloc, head, 1);
    }
};
