#include<iostream>
#include<vector>

template<typename T, typename Allocator = std::allocator<T>>
class List {
private:
    struct Node {
        T value;
        Node* prev;
        Node* next;

        explicit Node(const T& value, Node* prev, Node* next) :value(value), prev(prev), next(next) {};
        explicit Node(T&& value, Node* prev, Node* next) :value(std::move(value)), prev(prev), next(next) {};

        template<typename U, typename V>
        Node(U&& first, V&& second, Node* prev, Node* next) : value(std::move(first), std::move(second)), prev(prev),
                next(next) {}

        Node(Node* prev = nullptr, Node* next = nullptr): value(T()), prev(prev), next(next) {}
    };

    using NodeAlloc = typename std::allocator_traits<Allocator>::template rebind_alloc<Node>;
public:
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

        common_iterator(const Node* node) : node(const_cast<Node *>(node)) {}
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

    List(List&& other) {
        head = other.head;
        other.head = nullptr;
        sz = other.sz;
        other.sz = 0;
        node_alloc = other.node_alloc;
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

    List& operator=(List&& other) {
        if (this == &other) {
            other.clear();
            return *this;
        }
        if(head) {
            auto current_node = head->next;
            while (current_node != head) {
                auto next = current_node->next;
                destroy_deallocate_node(current_node);
                current_node = next;
            }
            head->prev = head;
            head->next = head;
        }
        sz = 0;
        head = other.head;
        other.head = nullptr;
        sz = other.sz;
        other.sz = 0;
        node_alloc = other.node_alloc;
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

    void insert_before(const_iterator pos, T&& value) {
        auto prev = std::prev(pos).node;
        auto new_node = std::allocator_traits<NodeAlloc>::allocate(node_alloc, 1);
        std::allocator_traits<NodeAlloc>::construct(node_alloc, new_node, std::move(value), prev, pos.node);
        link_nodes(new_node, pos.node);
        link_nodes(prev, new_node);
        ++sz;
    }

    template<typename U, typename V>
    void insert_before(const_iterator pos, U&& first, V&& second) {
        auto prev = std::prev(pos).node;
        auto new_node = std::allocator_traits<NodeAlloc>::allocate(node_alloc, 1);
        std::allocator_traits<NodeAlloc>::construct(node_alloc, new_node, std::move(first), std::move(second), prev, pos.node);
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

    void clear() {
        auto current_node = head->next;
        while (current_node != head) {
            auto next = current_node->next;
            destroy_deallocate_node(current_node);
            current_node = next;
        }
        head->prev = head->next = head;
        sz = 0;
    }

    ~List() {
        if(head) {
            auto current_node = head->next;
            while (current_node != head) {
                auto next = current_node->next;
                destroy_deallocate_node(current_node);
                current_node = next;
            }
            head->prev = head->next = head;
        }
        sz = 0;
        std::allocator_traits<NodeAlloc>::deallocate(node_alloc, head, 1);
    }
};



template<typename Key, typename Value, typename Hash = std::hash<Key>, typename Equal = std::equal_to<Key>, typename Alloc = std::allocator<std::pair<const Key, Value>>>
class UnorderedMap {
public:
    using NodeType = std::pair<const Key, Value>;

    //-------Iterators--------//

    using Iterator = typename List<NodeType, Alloc>::iterator;
    using ConstIterator = typename List<NodeType, Alloc>::const_iterator;

    Iterator begin() { return elements.begin(); }
    ConstIterator begin() const { return ConstIterator(elements.head->next); }
    ConstIterator cbegin() const { return ConstIterator(elements.head->next); }
    Iterator end() { return elements.end(); }
    ConstIterator end() const { return ConstIterator(elements.head); }
    ConstIterator cend() const { return ConstIterator(elements.head); }

    //---------------------//

private:

    List<NodeType, Alloc> elements;
    std::vector<List<Iterator, Alloc>> hash_table;
    size_t table_size = 811;
    size_t sz = 0;

    Hash hash =  Hash();
    Equal equal = Equal();
    Alloc alloc = Alloc();


public:

//  UnorderedMap() = default;

    UnorderedMap(size_t table_size = 811, const Hash& hash = Hash(), const Equal& equal = Equal(), const Alloc& alloc = Alloc()) :
            table_size(table_size), hash(hash), equal(equal), alloc(alloc) {
        hash_table.resize(table_size);
    }

    UnorderedMap(const UnorderedMap& other)  {
        //elements = List<NodeType, Alloc>(other.elements);
        elements = other.elements;
        sz = other.sz;
        hash_table = other.hash_table;
        table_size = other.table_size;
        hash = Hash();
        equal = Equal();
        alloc = Alloc();
    }

    UnorderedMap(UnorderedMap&& other) :
            elements(std::move(other.elements)), hash_table(std::move(other.hash_table)), table_size(other.table_size) ,
            sz(other.sz) {
        other.sz = 0;
    }

    UnorderedMap& operator=(const UnorderedMap& other) {
        elements = List<NodeType, Alloc>(other.elements);
        sz = other.sz;
        hash_table = other.hash_table;
        table_size = other.table_size;
        hash = Hash();
        equal = Equal();
        alloc = Alloc();
        return *this;
    }

    UnorderedMap& operator=(UnorderedMap&& other) {
        elements = std::move(other.elements);
        sz = other.sz;
        other.sz = 0;
        hash_table = std::move(other.hash_table);
        table_size = other.table_size;
        other.table_size = 811;
        hash = other.hash;
        equal = other.equal;
        alloc = other.alloc;
        return *this;
    }

    Value& operator[](Key key) {
        size_t h = hash(key) % table_size;
        for(auto& it : hash_table[h]) {
            if(equal((*it).first, key))
                return (*it).second;
        }
        if(load_factor() >= max_load_factor()) {
            reserve(6 * table_size + 1);
            h = hash(key) % table_size;
        }
        ++sz;
        elements.insert_before(elements.begin(), NodeType(key, Value()));
        hash_table[h].insert_before(hash_table[h].end(), elements.begin());
        return (*elements.begin()).second;

    }

    Value& at(Key key) {
        size_t h = hash(key) % table_size;
        for(auto& it : hash_table[h]) {
            if(equal((*it).first, key))
                return (*it).second;
        }
        throw std::out_of_range("Sorry, not in map");
    }

    const Value& at(const Key key) const {
        size_t h = hash(key) % table_size;
        for(auto& it : hash_table[h]) {
            if(equal((*it).first, key))
                return (*it).second;
        }
        throw std::out_of_range("Sorry, not in map");
    }

    size_t size() { return sz;}

    std::pair<Iterator, bool> insert(const NodeType& node) {
        size_t h = hash(node.first) % table_size;
        for(auto& it : hash_table[h]) {
            if(equal((*it).first, node.first))
                return std::make_pair(it, false);
        }
        if(load_factor() >= max_load_factor()) {
            reserve(6 * table_size + 1);
            h = hash(node.first) % table_size;
        }
        ++sz;
        elements.insert_before(elements.begin(), std::move(node));
        hash_table[h].insert_before(hash_table[h].end(), elements.begin());
        return std::make_pair(elements.begin(), true);
    }

    auto insert(NodeType&& node) -> std::pair<Iterator, bool> {
        size_t h = hash(node.first) % table_size;
        for(auto& it : hash_table[h]) {
            if(equal((*it).first, node.first))
                return std::make_pair(it, false);
        }
        if(load_factor() >= max_load_factor()) {
            reserve(6 * table_size + 1);
            h = hash(node.first) % table_size;
        }
        ++sz;
        elements.insert_before(elements.begin(), std::move(const_cast<Key&>(node.first)), std::move(node.second));
        hash_table[h].insert_before(hash_table[h].end(), elements.begin());
        return std::make_pair(elements.begin(), true);
    }

    template <typename InputIterator>
    void insert(InputIterator begin, InputIterator end) {
        for(auto it = begin; it != end; ++it) {
            insert(*it);
        }
    }

    template< class... Args >
    auto emplace(Args&&... args) -> std::pair<Iterator, bool> {
        auto new_node = std::allocator_traits<Alloc>::allocate(alloc, 1);
        std::allocator_traits<Alloc>::construct(alloc, new_node, std::forward<Args>(args)...);
        std::pair<Iterator, bool> res = insert(std::move(*new_node));
        std::allocator_traits<Alloc>::destroy(alloc, new_node);
        std::allocator_traits<Alloc>::deallocate(alloc, new_node, 1);
        return res;
    }

    void erase(Iterator it) {
        size_t h = hash((*it).first) % table_size;
        for(auto i = hash_table[h].begin(); i != hash_table[h].end(); ++i) {
            if((*i) == it)
                hash_table[h].erase(i);
        }
        --sz;
        elements.erase(it);
    }

    void erase(Iterator begin, Iterator end) {
        for(auto it = begin; it != end; ++it) {
            erase(it);
        }
    }

    Iterator find(Key key) {
        size_t h = hash(key) % table_size;
        for(auto& it : hash_table[h]) {
            if(equal((*it).first, key))
                return it;
        }
        return elements.end();
    }

    void reserve(size_t count) {
        std::vector<List<Iterator, Alloc>> new_hash_table;
        new_hash_table.resize(count);
        for(const auto& row : hash_table) {
            for(const auto& elem : row) {
                size_t h = hash((*elem).first) % count;
                new_hash_table[h].insert_before(new_hash_table[h].end(), elem);
            }
        }
        hash_table = std::move(new_hash_table);
        table_size = count;
    }

    size_t max_size() const noexcept {
        return table_size;
    }

    float load_factor() const {
        return static_cast<float>(sz) / static_cast<float>(table_size);
    }

    float max_load_factor() const {
        return 1.0;
    }

    ~UnorderedMap() = default;

    friend class List<NodeType, Alloc>;

};
