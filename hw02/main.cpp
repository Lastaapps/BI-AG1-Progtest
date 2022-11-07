#ifndef __PROGTEST__
#include <cassert>
#include <iostream>
#include <memory>
#include <limits>
#include <optional>
#include <algorithm>
#include <bitset>
#include <list>
#include <array>
#include <vector>
#include <deque>
#include <unordered_set>
#include <unordered_map>
#include <stack>
#include <queue>
#include <random>

#endif

using namespace std;

template<typename T>
ostream& operator<<(ostream& out, const shared_ptr<const T>& t) noexcept {
  return (t == nullptr) ? out << "null" : out << *t;
}

template<typename U>
size_t safeSum(const shared_ptr<const U>& node) noexcept {
  return (node == nullptr) ? 0 : node -> sum;
}

template<typename U>
size_t safeProdSum(const shared_ptr<const U>& node) noexcept {
  return (node == nullptr) ? 0 : node -> sellSum;
}

template<typename Key>
struct BaseNode {
  using SharedNode = shared_ptr<BaseNode>;

  Key key;
  size_t sum = 1;
  char weight = 0;
  SharedNode left = nullptr, right = nullptr;
  weak_ptr<BaseNode> parent;

  BaseNode(const Key& k) noexcept : key(k) {};

  bool isLeftChild() const noexcept {
    const SharedNode lock = parent.lock();
    return lock == nullptr || lock -> left.get() == this;
  }

  bool isList() const noexcept {
    return left == nullptr && right == nullptr;
  }

  bool isFull() const noexcept {
    return left != nullptr && right != nullptr;
  }

  void reSum() noexcept {
    sum = safeSum<BaseNode>(left) + safeSum<BaseNode>(right) + 1;
  }

  friend ostream& operator<<(ostream& out, const BaseNode& node) noexcept {
    const SharedNode p = node.parent.lock();
    out << "{"
      << node.key
      << ", w: " << (int)node.weight
      << ", s: " << node.sum
      << ", l: " << node.isLeftChild()
      << ", p: ";
    (p == nullptr ? out << "null"s : out << p -> key)
      << "}";
    return out;
  }

  friend void swapKey(BaseNode& n0, BaseNode& n1) noexcept { 
    swap(n0.key, n1.key);
    n0.reSum();
    n1.reSum();
  }
};


template <typename Prod>
struct Holder {
  const size_t items = 0;
  const Prod data;

  bool operator<(const Holder& other) const noexcept {
    return items < other.items ? true 
      : items > other.items ? false 
      : data < other.data;
  }
};

template <typename Key>
struct ProdNode final : public BaseNode<Holder<Key>> {
  using HKey = Holder<Key>;
  size_t sellSum;

  ProdNode(const Key& key) noexcept : BaseNode<HKey>(key), sellSum(key.items) {}

  void reSum() noexcept {
    BaseNode<HKey>::reSum();
    sellSum = BaseNode<HKey>::key.items + safeProdSum<HKey>(left) + safeProdSum<HKey>(right);
  }
};

template <typename Key, typename Node = BaseNode<Key>, typename Cmp = less<Key>>
class AVLTree {

  using SharedNode = shared_ptr<Node>;
  SharedNode root = nullptr;

  struct Iterator {
    SharedNode ptr;
    Node* operator->() noexcept {
      return ptr.get();
    }
    Node& operator*() noexcept {
      return *ptr;
    }
    void operator++() noexcept {
      if (ptr -> right == nullptr) {
        while (ptr != nullptr) {
          const bool isNextLeft = ptr -> isLeftChild();
          ptr = ptr -> parent.lock();
          if (isNextLeft) break;
        }
      } else {
        ptr = ptr -> right;
        while (true) {
          SharedNode next = ptr -> left;
          if (next == nullptr) break;
          ptr = next;
        }
      }
    }
    bool operator==(const Iterator& other) noexcept {
      return ptr == other.ptr;
    }
    bool operator!=(const Iterator& other) noexcept {
      return ptr != other.ptr;
    }
  };

  public:
  Iterator begin() noexcept {
    if (root == nullptr)
      return Iterator({nullptr});

    SharedNode node = root;
    while (true) {
      SharedNode next = node -> left;
      if (next -> left == nullptr)
        break;
      node = next;
    }
    return Iterator({node});
  }

  Iterator end() const noexcept {
    return Iterator({nullptr});
  }

  Iterator byIndexImpl(const size_t k, SharedNode node) noexcept {
    if (node == nullptr) return end();

    const size_t lSum = safeSum<Node>(node -> left);
    if (k < lSum) {
      return byIndexImpl(k, node -> left);
    } else if (k == lSum) {
      return Iterator{node};
    } else {
      return byIndexImpl(k - lSum - 1, node -> right);
    }
  }

  Iterator byIndex(const size_t index) noexcept {
    return byIndexImpl(index, root);
  }

  size_t indexOfImpl(SharedNode node) {
    SharedNode parent = node -> parent.lock();
    if (parent == nullptr) return 0;

    const bool isLeft = node -> isLeftChild();
    return (isLeft ? 0 : safeSum<Node>(parent) - safeSum<Node>(node)) + indexOfImpl(parent);
  }

  size_t indexOf(Iterator iter) noexcept {
    SharedNode node = iter.ptr;
    return safeSum<Node>(node -> left) + indexOfImpl(node);
  }

  private:
  void setChild(SharedNode parent, SharedNode child, bool isLeft) noexcept {
    if (parent != nullptr)
      (isLeft ? parent -> left : parent -> right) = child;
    else
      root = child;

    if (child != nullptr)
      child -> parent = parent;
  }

  SharedNode rotateNode(SharedNode node) noexcept {
    const char weight = node -> weight;
    SharedNode parent = node -> parent.lock();
    const bool isLeft = node -> isLeftChild();

    if (weight == 2) {
      SharedNode n1 = node;
      SharedNode n2 = n1 -> right;
      // right rotation
      if (n2 -> weight != -1) {
        // simple rotation
        cout << "Simple right rotation" << endl;
        SharedNode a  = n1 -> left;
        SharedNode b  = n2 -> left;
        SharedNode c  = n2 -> right;

        if (n2 -> weight == 0) {
          n1 -> weight = 1;
          n2 -> weight = -1;
        } else { // 1
          n1 -> weight = 0;
          n2 -> weight = 0;
        }

        setChild(parent, n2, isLeft);
        setChild(n2, n1, true);
        setChild(n2,  c, false);
        setChild(n1,  a, true);
        setChild(n1,  b, false);

        n1 -> reSum();
        n2 -> reSum();

        return n2;

      } else {
        // double rotation
        cout << "Double right rotation" << endl;
        SharedNode n1 = node;
        SharedNode n2 = node -> right;
        SharedNode n3 = n2 -> left;
        SharedNode a  = n1 -> left;
        SharedNode b  = n3 -> left;
        SharedNode c  = n3 -> right;
        SharedNode d  = n2 -> right;

        n1 -> weight = n3 -> weight == 1 ? 1 : 0;
        n2 -> weight = n3 -> weight == -1 ? -1 : 0;
        n3 -> weight = 0;

        setChild(parent, n3, isLeft);
        setChild(n3, n1, true);
        setChild(n3, n2, false);
        setChild(n1,  a, true);
        setChild(n1,  b, false);
        setChild(n2,  c, true);
        setChild(n2,  d, false);

        n1 -> reSum();
        n2 -> reSum();
        n3 -> reSum();

        return n3;
      }
    } else if (weight == -2) {
      SharedNode n1 = node;
      SharedNode n2 = n1 -> left;
      if (n2 -> weight != 1) {
        // simple rotation
        cout << "Simple left rotation" << endl;
        SharedNode a  = n1 -> right;
        SharedNode b  = n2 -> right;
        SharedNode c  = n2 -> left;

        if (n2 -> weight == 0) {
          n1 -> weight = -1;
          n2 -> weight = 1;
        } else { // 1
          n1 -> weight = 0;
          n2 -> weight = 0;
        }

        setChild(parent, n2, isLeft);
        setChild(n2, n1, false);
        setChild(n2,  c, true);
        setChild(n1,  a, false);
        setChild(n1,  b, true);

        n1 -> reSum();
        n2 -> reSum();

        return n2;

      } else {
        // double rotation
        cout << "Double left rotation" << endl;
        SharedNode n1 = node;
        SharedNode n2 = node -> left;
        SharedNode n3 = n2 -> right;
        SharedNode a  = n1 -> right;
        SharedNode b  = n3 -> right;
        SharedNode c  = n3 -> left;
        SharedNode d  = n2 -> left;

        n1 -> weight = n3 -> weight == -1 ? -1 : 0;
        n2 -> weight = n3 -> weight == 1 ? 1 : 0;
        n3 -> weight = 0;

        setChild(parent, n3, isLeft);
        setChild(n3, n1, false);
        setChild(n3, n2, true);
        setChild(n1,  a, false);
        setChild(n1,  b, true);
        setChild(n2,  c, false);
        setChild(n2,  d, true);

        n1 -> reSum();
        n2 -> reSum();
        n3 -> reSum();

        return n3;
      }
    } else {
#ifndef __PROGTEST__
      throw runtime_error("You fucked up!");
#endif
    }
  }

  bool insertImpl(const Key& key, SharedNode node) noexcept {

    bool isLeft = Cmp{}(key, node -> key);
    SharedNode child = isLeft ? node -> left : node -> right;

    if (child == nullptr) {
      SharedNode toInsert = make_shared<Node>(key);
      setChild(node, toInsert, isLeft);
      cout << "Incing   " << *node << endl;
      node -> weight += isLeft ? -1 : 1;
      node -> reSum();
      return node -> weight != 0;
    } else {
      const bool incLevel = insertImpl(key, child);
      node -> reSum();
      if (incLevel) {
        node -> weight += isLeft ? -1 : 1;
        cout << "Incing   " << *node << endl;
        if (abs(node -> weight) == 2) {
          cout << "Rotating " << *node << endl;
          node = rotateNode(node);
          cout << "Rot res  " << *node << endl;
        }
      }
      return node -> weight != 0 && incLevel;
    }
  }

  public:
  void insert(const Key& key) noexcept {
    cout << "Insert   " << key << endl;
    SharedNode n = root;

    if (n == nullptr) {
      root = make_shared<Node>(key);
      return;
    }

    insertImpl(key, n);
  }
  private:
  void sendDecreaseUp(SharedNode node, bool fromLeft, bool sumOnly) noexcept {
    if (node == nullptr) return;

    SharedNode parent = node -> parent.lock();
    node -> reSum();

    if (sumOnly)
      return sendDecreaseUp(parent, false, true);

    node -> weight += -1 * (fromLeft ? -1 : 1);

    cout << "Decrease " << *node << endl;
    if (abs(node -> weight) == 2) {
      cout << "Rotating " << *node << endl;
      node = rotateNode(node);
      cout << "Rot res  " << *node << endl;

      if (node -> weight == 0)
        sendDecreaseUp(parent, node -> isLeftChild(), false);
    } 
    else if (node -> weight == 0)
      sendDecreaseUp(parent, node -> isLeftChild(), false);
    else
      sendDecreaseUp(parent, false, true);
  }

  public:
  void deleteNode(Iterator& iter) noexcept {
    SharedNode node = iter.ptr;
    cout << "Deleting " << iter -> key << endl;
    bool isLeft = node -> isLeftChild();
    SharedNode parent = node -> parent.lock();

    if (node -> isList()) {
      bool isLeft = node -> isLeftChild();
      SharedNode parent = node -> parent.lock();

      setChild(parent, SharedNode(), isLeft);
      sendDecreaseUp(parent, isLeft, false);

    } else if (node -> isFull()) {
      Iterator sucIter = iter; ++sucIter;
      SharedNode suc = sucIter.ptr;
      swapKey(*node, *suc);
      deleteNode(sucIter);
    } else {
      SharedNode child = node -> left != nullptr ? node -> left : node -> right;

      setChild(parent, child, isLeft);
      sendDecreaseUp(parent, isLeft, false);
    }
  }

  private:
  Iterator lowerBoundImpl(const Key& key, const SharedNode node) const noexcept {
    if (node == nullptr) return end();

    cout << "Lowering " << *node << endl;
    if (!Cmp{}(node -> key, key)) {
      if (!Cmp{}(key, node -> key)) return Iterator{node};

      const SharedNode left = node -> left;
      if (left == nullptr) return Iterator{node};

      return lowerBoundImpl(key, left);
    } else {
      return lowerBoundImpl(key, node -> right);
    }
  }

  public:
  Iterator lowerBound(const Key& key) const noexcept {
    cout << "Lower    " << key << endl;
    return lowerBoundImpl(key, root);
  }

  private:
  Iterator upperBoundImpl(const Key& key, const SharedNode node) const noexcept {
    if (node == nullptr) return end();

    cout << "Uppering " << *node << endl;
    if (!Cmp{}(key, node -> key)) {
      if (!Cmp{}(node -> key, key)) return Iterator{node};

      const SharedNode right = node -> right;
      if (right == nullptr) return Iterator{node};

      return upperBoundImpl(key, right);
    } else {
      return upperBoundImpl(key, node -> left);
    }
  }

  public:
  Iterator upperBound(const Key& key) const noexcept {
    cout << "Upper    " << key << endl;
    return upperBoundImpl(key, root);
  }

  size_t size() const noexcept {
    return root != nullptr ? root -> sum : 0;
  }

  void printImpl(ostream& out, SharedNode node, size_t depth) const noexcept {
    if (node == nullptr) return;
    for (size_t i = 0;  i < depth; ++i) out << " ";
    out << *node << endl;
    printImpl(out, node -> left, depth + 1);
    printImpl(out, node -> right, depth + 1);
  }

  void print(ostream& out = cout) const noexcept { printImpl(out, root, 0); }
};

template<typename Product>
struct Bestsellers {
  using PHolder = Holder<Product>;
  AVLTree<PHolder, ProdNode<Product>, greater<PHolder>> tree;

  // The total number of tracked products
  size_t products() const;

  void sell(const Product& p, size_t amount);

  // The most sold product has rank 1
  size_t rank(const Product& p) const;
  const Product& product(size_t rank) const;

  // How many copies of product with given rank were sold
  size_t sold(size_t rank) const;
  // The same but sum over interval of products (including from and to)
  // It must hold: sold(x) == sold(x, x)
  size_t sold(size_t from, size_t to) const;

  // Bonus only, ignore if you are not interested in bonus
  // The smallest (resp. largest) rank with sold(rank) == sold(r)
  size_t first_same(size_t r) const { return 0; }
  size_t last_same(size_t r) const { return 0; }
};

#ifndef __PROGTEST__

void treeTest() {
  // Right simple rotation
  //{
  //  AVLTree<ssize_t> tree;
  //  for (size_t i = 0; i < 20; ++i) {
  //    tree.insert(i);
  //    tree.print();
  //  }
  //  cout << "\n\n\n";
  //}

  // Left simple rotation
  //{
  //  AVLTree<ssize_t> tree;
  //  for (size_t i = 20; i > 0; --i) {
  //    tree.insert(i - 1);
  //    tree.print();
  //  }
  //  cout << "\n\n\n";
  //}
  //{

  // Double right/left rotation
  //{
  //  AVLTree<ssize_t> tree;
  //  for (size_t i = 0; i < 20; ++i) {
  //    tree.insert(i % 2 == 0 ? i / 2 : 19 - i / 2);
  //    tree.print();
  //  }
  //  cout << "\n\n\n";
  //}

  // Iterator
  //  AVLTree<ssize_t> tree;
  //  for (size_t i = 20; i > 0; --i) tree.insert(i - 1);
  //  tree.print();
  //   for (auto itr = tree.begin(); itr != tree.end(); ++itr)
  //     cout << *itr << endl;
  //  cout << "\n\n\n";
  //}

  // Delete
  //{
  //  AVLTree<ssize_t> tree;
  //  for (size_t i = 0; i < 20; ++i)
  //    tree.insert(i);
  //  tree.print();
  //  for (size_t i = 1; i < 20; i += 3) {
  //    auto bound = tree.lowerBound(i);
  //    tree.deleteNode(bound);
  //    tree.print();
  //  }
  //  cout << "\n\n\n";
  //}

  // By index
  {
    AVLTree<ssize_t> tree;
    for (size_t i = 0; i < 20; ++i)
      tree.insert(i);
    tree.print();
    for (size_t i = 0; i < 20; ++i)
      cout << i << ": " << *tree.byIndex(i) << endl;
    for (size_t i = 0; i < 20; ++i)
      cout << i << ": " << tree.indexOf(tree.byIndex(i)) << endl;
    cout << "\n\n\n";
  }
}

void test1() {
  Bestsellers<std::string> T;
//   T.sell("coke", 32);
//   T.sell("bread", 1);
//   assert(T.products() == 2);
//   T.sell("ham", 2);
//   T.sell("mushrooms", 12);
// 
//   assert(T.products() == 4);
//   assert(T.rank("ham") == 3);
//   assert(T.rank("coke") == 1);
//   assert(T.sold(1, 3) == 46);
//   assert(T.product(2) == "mushrooms");
// 
//   T.sell("ham", 11);
//   assert(T.products() == 4);
//   assert(T.product(2) == "ham");
//   assert(T.sold(2) == 13);
//   assert(T.sold(2, 2) == 13);
//   assert(T.sold(1, 2) == 45);
}
 
// void test2() {
// # define CATCH(expr) \
//   try { expr; assert(0); } catch (const std::out_of_range&) { assert(1); };
// 
//   Bestsellers<std::string> T;
//   T.sell("coke", 32);
//   T.sell("bread", 1);
// 
//   CATCH(T.rank("ham"));
//   CATCH(T.product(3));
//   CATCH(T.sold(0));
//   CATCH(T.sold(9));
//   CATCH(T.sold(0, 1));
//   CATCH(T.sold(3, 2));
//   CATCH(T.sold(1, 9));
// 
// #undef CATCH
// }

int main() {
  treeTest();
  test1();
  // test2();
  return 0;
}

#endif

