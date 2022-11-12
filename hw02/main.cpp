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

#define COUT if constexpr (true) cout

#elif

#define COUT if constexpr (false) cout

#endif

using namespace std;

template<typename T>
ostream& operator<<(ostream& out, const shared_ptr<const T>& t) noexcept {
  return (t == nullptr) ? out << "null" : out << *t;
}


// --- BaseNode ---------------------------------------------------------------
/** Gets nodes sum or 0 for an empty node*/
template<typename U>
size_t safeSum(const shared_ptr<const U>& node) noexcept {
  return (node == nullptr) ? 0 : node -> sum;
}

/** Base node implementation so I can keep  this for the future */
template<typename Key>
struct BaseNode {
  using SharedNode = shared_ptr<BaseNode>;

  Key key;
  // number of children plus 1 for this node
  size_t sum = 1;
  // +1 if right node has more levels, -1 for left and 0 for equality
  char weight = 0;
  // children in a tree
  SharedNode left = nullptr, right = nullptr;
  // a parent in a tree
  weak_ptr<BaseNode> parent;

  BaseNode(const Key& k) noexcept : key(k) {};
  virtual ~BaseNode(){};

  /** If this node is a left child of it's parent of it is root */
  bool isLeftChild() const noexcept {
    const SharedNode lock = parent.lock();
    return lock == nullptr || lock -> left.get() == this;
  }

  /** If the node is list, e.g. has no children */
  bool isList() const noexcept {
    return left == nullptr && right == nullptr;
  }

  /** If the node has 2 children */
  bool isFull() const noexcept {
    return left != nullptr && right != nullptr;
  }

  /** Recomputes sum property */
  virtual void reSum() noexcept {
    sum = safeSum<BaseNode>(left) + safeSum<BaseNode>(right) + 1;
  }

  /** Prints the node out */
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

  // product count may have to be recalculated, not in the base class
  virtual void onKeyChanged() noexcept {}

  /** Exchanges key with another node */
  friend void swapKey(BaseNode& n0, BaseNode& n1) noexcept { 
    swap(n0.key, n1.key);
    n0.onKeyChanged();
    n1.onKeyChanged();
  }
};



// --- Holder -----------------------------------------------------------------
/** Compares nodes according to this assignment, using greatest */
template <typename Holder>
struct HolderCmp final {
  bool operator()(const Holder& h0, const Holder& h1) const noexcept {
    return
      h0.items > h1.items ? true :
      h0.items < h1.items ? false :
      h0.isLowest < h1.isLowest ? true :
      h0.isLowest > h1.isLowest ? false :
      h0.data < h1.data;
  }
};

/** Gets nodes product count or 0 for an empty node */
template<typename U>
size_t safeProdSum(const shared_ptr<const U>& node) noexcept {
  return (node == nullptr) ? 0 : node -> sellSum;
}

/** Node content for this assignment */
template <typename Prod>
struct Holder final {
  // number of items sold
  const size_t items = 0;
  // to make lower/upper bound easier
  const bool isLowest = false;
  // sold product
  const Prod data;
};

/** Holds holder in the tree, adds sold sums */
template <typename Prod>
struct ProdNode final : public BaseNode<Holder<Prod>> {
  using HKey = Holder<Prod>;
  // number of items sold in subtrees and in this node
  size_t sellSum;

  ProdNode(const Prod& key) noexcept : BaseNode<HKey>(key), sellSum(key.items) {}

  /** Also resums the number of products sold */
  void reSum() noexcept override {
    BaseNode<HKey>::reSum();
    sellSum = BaseNode<HKey>::key.items + safeProdSum<HKey>(left) + safeProdSum<HKey>(right);
  }

  void onKeyChanged() noexcept override {
    // product count may have to be recalculated
    reSum();
    // update patents
    const auto p = BaseNode<HKey>::parent.lock();
    if (p != nullptr) p -> reSum();
  }
};



// --- AVL tree ---------------------------------------------------------------
template <typename Key, typename Node = BaseNode<Key>, typename Cmp = less<Key>>
class AVLTree {

  using SharedNode = shared_ptr<Node>;
  using SharedConstNode = shared_ptr<const Node>;

  // tree root node
  SharedNode root = nullptr;

  // comparator instance
  const Cmp cmp = {};


  // iterates ower the tree
  struct Iterator {
    SharedConstNode ptr;

    SharedNode unsafePtr() noexcept {
      return const_pointer_cast<Node>(ptr);
    }
    Node* operator->() noexcept { return unsafePtr().get(); }
    Node& operator*() noexcept { return *unsafePtr(); }

    const Node* operator->() const noexcept { return ptr.get(); }
    const Node& operator*() const noexcept { return *ptr; }

    bool operator==(const Iterator& other) const noexcept { return ptr == other.ptr; }
    bool operator!=(const Iterator& other) const noexcept { return ptr != other.ptr; }

    /** successor */
    void operator++() noexcept {
      if (ptr -> right == nullptr) {
        // find the first parent in the right
        while (ptr != nullptr) {
          const bool isNextLeft = ptr -> isLeftChild();
          ptr = ptr -> parent.lock();
          if (isNextLeft) break;
        }
      } else {
        // find the first on the right and then go left
        ptr = ptr -> right;
        while (true) {
          SharedNode next = ptr -> left;
          if (next == nullptr) break;
          ptr = next;
        }
      }
    }
  };

  public:
  /** Iterator to the first tree item */
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

  /** Iterator after the end of the tree */
  Iterator end() const noexcept { return Iterator({nullptr}); }



  private:
  /** Gets item in the tree by it's index */
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

  public:
  /** Gets item in the tree by it's index */
  Iterator byIndex(const size_t index) noexcept {
    return byIndexImpl(index, root);
  }



  private:
  /** Gets index of an item in the tree */
  size_t indexOfImpl(SharedNode node) {
    SharedNode parent = node -> parent.lock();
    if (parent == nullptr) return 0;

    const bool isLeft = node -> isLeftChild();
    return (isLeft ? 0 : safeSum<Node>(parent) - safeSum<Node>(node)) + indexOfImpl(parent);
  }

  public:
  /** Gets index of an item in the tree */
  size_t indexOf(Iterator iter) noexcept {
    SharedNode node = iter.ptr;
    return safeSum<Node>(node -> left) + indexOfImpl(node);
  }



  private:
  /** Connects parent and a child together.
   * No recalculation is done.
   * If parent is nullptr, child is set as the tree root
   */
  void setChild(SharedNode parent, SharedNode child, bool isLeft) noexcept {
    if (parent != nullptr)
      (isLeft ? parent -> left : parent -> right) = child;
    else
      root = child;

    if (child != nullptr)
      child -> parent = parent;
  }

  /** Handles all the awl rotations based on the node weight */
  SharedNode rotateNode(SharedNode node) noexcept {

    const char weight = node -> weight;
    const SharedNode parent = node -> parent.lock();
    const bool isLeft = node -> isLeftChild();

    if (weight == 2) {

      const SharedNode n1 = node;
      const SharedNode n2 = n1 -> right;

      // right rotations
      if (n2 -> weight != -1) {
        // simple rotation
        COUT << "Simple right rotation" << endl;
        const SharedNode a  = n1 -> left;
        const SharedNode b  = n2 -> left;
        const SharedNode c  = n2 -> right;

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
        COUT << "Double right rotation" << endl;
        const SharedNode n3 = n2 -> left;
        const SharedNode a  = n1 -> left;
        const SharedNode b  = n3 -> left;
        const SharedNode c  = n3 -> right;
        const SharedNode d  = n2 -> right;

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

    // left rotations
    } else if (weight == -2) {

      const SharedNode n1 = node;
      const SharedNode n2 = n1 -> left;

      if (n2 -> weight != 1) {
        // simple rotation
        COUT << "Simple left rotation" << endl;
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
        COUT << "Double left rotation" << endl;
        const SharedNode n3 = n2 -> right;
        const SharedNode a  = n1 -> right;
        const SharedNode b  = n3 -> right;
        const SharedNode c  = n3 -> left;
        const SharedNode d  = n2 -> left;

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



  private:
  /** Inserts a new item into the tree, key mas be unique */
  bool insertImpl(const Key& key, SharedNode node) noexcept {

    const bool isLeft = cmp(key, node -> key);
    const SharedNode child = isLeft ? node -> left : node -> right;

    if (child == nullptr) {
      SharedNode toInsert = make_shared<Node>(key);
      setChild(node, toInsert, isLeft);
      COUT << "Incing   " << *node << endl;
      node -> weight += isLeft ? -1 : 1;
      node -> reSum();
      return node -> weight != 0;
    } else {
      const bool incLevel = insertImpl(key, child);
      node -> reSum();
      if (incLevel) {
        node -> weight += isLeft ? -1 : 1;
        COUT << "Incing   " << *node << endl;
        if (abs(node -> weight) == 2) {
          COUT << "Rotating " << *node << endl;
          node = rotateNode(node);
          COUT << "Rot res  " << *node << endl;
        }
      }
      return node -> weight != 0 && incLevel;
    }
  }

  public:
  /** Inserts a new item into the tree, key mas be unique */
  void insert(const Key& key) noexcept {
    COUT << "Insert   " << key << endl;
    const SharedNode n = root;

    if (n == nullptr)
      // insert root
      root = make_shared<Node>(key);
    else
      // insert new node
      insertImpl(key, n);
  }



  private:
  /** Propagates info about node deletion to it's parents */
  void sendDecreaseUp(const SharedNode node, const bool fromLeft, const bool sumOnly) noexcept {
    if (node == nullptr) return;

    const SharedNode parent = node -> parent.lock();
    node -> reSum();

    if (sumOnly)
      return sendDecreaseUp(parent, false, true);

    // add node decrease info, inverse of insert
    node -> weight += -1 * (fromLeft ? -1 : 1);

    COUT << "Decrease " << *node << endl;
    if (abs(node -> weight) == 2) {
      COUT << "Rotating " << *node << endl;
      node = rotateNode(node);
      COUT << "Rot res  " << *node << endl;

      if (node -> weight == 0)
        sendDecreaseUp(parent, node -> isLeftChild(), false);
    } 
    else if (node -> weight == 0)
      sendDecreaseUp(parent, node -> isLeftChild(), false);
    else
      // no further weights updates, just resuming
      sendDecreaseUp(parent, false, true);
  }



  public:
  // 
  void deleteNode(Iterator& iter) noexcept {
    const SharedNode node = iter.ptr;
    const bool isLeft = node -> isLeftChild();
    const SharedNode parent = node -> parent.lock();
    COUT << "Deleting " << iter -> key << endl;

    if (node -> isList()) {
      // Delete list
      const bool isLeft = node -> isLeftChild();
      const SharedNode parent = node -> parent.lock();

      setChild(parent, SharedNode(), isLeft);
      sendDecreaseUp(parent, isLeft, false);

    } else if (node -> isFull()) {
      // Delete node with 2 children
      Iterator sucIter = iter; ++sucIter;
      const SharedNode suc = sucIter.ptr;
      swapKey(*node, *suc);
      deleteNode(sucIter);
    } else {
      // Delete node with 1 child
      const SharedNode child = node -> left != nullptr ? node -> left : node -> right;

      setChild(parent, child, isLeft);
      sendDecreaseUp(parent, isLeft, false);
    }
  }
  


  private:
  /** Finds the lower bound for the key given */
  Iterator lowerBoundImpl(const Key& key, const SharedNode node) const noexcept {
    if (node == nullptr) return end();

    COUT << "Lowering " << *node << endl;
    if (!cmp(node -> key, key)) {
      if (!cmp(key, node -> key)) return Iterator{node};

      const SharedNode left = node -> left;
      if (left == nullptr) return Iterator{node};

      return lowerBoundImpl(key, left);
    } else {
      return lowerBoundImpl(key, node -> right);
    }
  }

  public:
  /** Finds the lower bound for the key given */
  Iterator lowerBound(const Key& key) const noexcept {
    COUT << "Lower    " << key << endl;
    return lowerBoundImpl(key, root);
  }



  private:
  /** Finds the item before upper bound */
  // TODO reimplement to support correct upper bound
  Iterator upperBoundImpl(const Key& key, const SharedNode node) const noexcept {
    if (node == nullptr) return end();

    COUT << "Uppering " << *node << endl;
    if (!cmp(key, node -> key)) {
      if (!cmp(node -> key, key)) return Iterator{node};

      const SharedNode right = node -> right;
      if (right == nullptr) return Iterator{node};

      return upperBoundImpl(key, right);
    } else {
      return upperBoundImpl(key, node -> left);
    }
  }

  public:
  /** Finds the item before upper bound */
  Iterator upperBound(const Key& key) const noexcept {
    COUT << "Upper    " << key << endl;
    return upperBoundImpl(key, root);
  }



  private:
  /** Finds the first common ancestor of two tree nodes */
  Iterator commonAncestorImpl(const Key& k0, const Key& k1, const SharedConstNode node) const noexcept {
    const Key& key = node -> key;
    const bool cmp0 = cmp(k0, key);
    const bool cmp1 = cmp(k1, key);
    COUT << "Ancestor " << key << " " << cmp0 << ":" << cmp1 << endl;

    if (cmp0 != cmp1 || k0 == key || k1 == key)
      return Iterator{node};

    return commonAncestorImpl(k0, k1, cmp0 ? node -> left : node -> right);
  }

  public:
  /** Finds the first common ancestor of two tree nodes */
  Iterator commonAncestor(const Iterator i0, const Iterator i1) const noexcept {
    COUT << "Common   " << i0.ptr << " " << i1.ptr << endl;
    if (i0 == end() || i1 == end()) return end();
    if (i0 == i1) return i0;
    return commonAncestorImpl(i0 -> key, i1 -> key, root);
  }


  /** Returns the number of nodes in the tree */
  size_t size() const noexcept { return root != nullptr ? root -> sum : 0; }



  private:
  /** Prints all the nodes in the tree with indentation representing the internal structure */
  void printImpl(ostream& out, SharedNode node, size_t depth) const noexcept {
    if (node == nullptr) return;
    for (size_t i = 0;  i < depth; ++i) out << " ";
    out << *node << endl;
    printImpl(out, node -> left, depth + 1);
    printImpl(out, node -> right, depth + 1);
  }

  public:
  /** Prints all the nodes in the tree with indentation representing the internal structure */
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
  //{
  //  AVLTree<ssize_t> tree;
  //  for (size_t i = 0; i < 20; ++i)
  //    tree.insert(i);
  //  tree.print();
  //  for (size_t i = 0; i < 20; ++i)
  //    cout << i << ": " << *tree.byIndex(i) << endl;
  //  for (size_t i = 0; i < 20; ++i)
  //    cout << i << ": " << tree.indexOf(tree.byIndex(i)) << endl;
  //  cout << "\n\n\n";
  //}

  // Common ancestor
  {
    AVLTree<ssize_t> tree;
    for (size_t i = 0; i < 20; ++i) tree.insert(i);
    tree.print();
    assert(tree.commonAncestor(tree.lowerBound(19), tree.lowerBound(6)) -> key == 7);
    assert(tree.commonAncestor(tree.lowerBound(14), tree.lowerBound(8)) -> key == 11);
    assert(tree.commonAncestor(tree.lowerBound(4), tree.lowerBound(6)) -> key == 5);
    assert(tree.commonAncestor(tree.lowerBound(19), tree.lowerBound(12)) -> key == 15);
    assert(tree.commonAncestor(tree.lowerBound(19), tree.lowerBound(17)) -> key == 17);
    assert(tree.commonAncestor(tree.lowerBound(17), tree.lowerBound(19)) -> key == 17);
    assert(tree.commonAncestor(tree.lowerBound(8), tree.lowerBound(8)) -> key == 8);
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

