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

#else

#define COUT if constexpr (false) cout

#endif

using namespace std;

// template<typename T>
// ostream& operator<<(ostream& out, const shared_ptr<const T>& t) noexcept {
//   return (t == nullptr) ? out << "null" : out << *t;
// }
// --- BaseNode ---------------------------------------------------------------
/** Gets nodes sum or 0 for an empty node*/
template<typename U>
size_t safeSum(const shared_ptr<const U>& node) noexcept {
  return (node == nullptr) ? 0 : node -> sum;
}

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
  const int8_t isLowest = 0;
  // sold product
  const Prod data;

  Holder(const size_t& i, const Prod& d)
    : items(i), data(d) {}
  Holder(const size_t& i, const int8_t l, const Prod& d)
    : items(i), isLowest(l), data(d) {}

  bool operator==(const Holder& other) const noexcept {
    return items == other.items && data == other.data;
  }
};

template <typename Prod>
ostream& operator<<(ostream& out, const Holder<Prod>& data) {
  return out << "<" << data.data << ", " << data.items << ">";
}

/** Base node implementation so I can keep  this for the future */
template<typename Key>
struct BaseNode final {
  using HKey = Holder<Key>;
  using SharedNode = shared_ptr<BaseNode>;

  const HKey key;
  // number of children plus 1 for this node
  size_t sum = 1;
  // +1 if right node has more levels, -1 for left and 0 for equality
  char weight = 0;
  // children in a tree
  SharedNode left = nullptr, right = nullptr;
  // a parent in a tree
  weak_ptr<BaseNode> parent;

  // redundant
  size_t sellSum;

  BaseNode(const HKey& k) noexcept : key(k), sellSum(k.items) {};
  ~BaseNode(){};

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
  void reSum() noexcept {
    sum = safeSum<BaseNode>(left) + safeSum<BaseNode>(right) + 1;
    sellSum = key.items + safeProdSum<BaseNode>(left) + safeProdSum<BaseNode>(right);
  }

  /** Prints the node out */
  friend ostream& operator<<(ostream& out, const BaseNode& node) noexcept {
    const SharedNode p = node.parent.lock();
    out << "{"
      << node.key
      << ", w: " << (int)node.weight
      << ", s: " << node.sum
      << ", d: " << node.sellSum
      << ", l: " << node.isLeftChild()
      << ", p: ";
    (p == nullptr ? out << "null"s : out << p -> key)
      << "}";
    return out;
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
  struct Iterator final {
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
  Iterator begin() const noexcept {
    if (root == nullptr)
      return Iterator({nullptr});

    SharedNode node = root;
    while (true) {
      SharedNode next = node -> left;
      if (next == nullptr)
        break;
      node = next;
    }
    return Iterator({node});
  }

  /** Iterator after the end of the tree */
  Iterator end() const noexcept { return Iterator({nullptr}); }



  private:
  /** Gets item in the tree by it's index */
  Iterator byIndexImpl(const size_t k, SharedConstNode node) const noexcept {
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
  Iterator byIndex(const size_t index) const noexcept {
    return byIndexImpl(index, root);
  }



  private:
  /** Gets index of an item in the tree */
  size_t indexOfImpl(const SharedConstNode node) const {
    SharedConstNode parent = node -> parent.lock();
    if (parent == nullptr) return 0;

    const bool isLeft = node -> isLeftChild();
    return (isLeft ? 0 : safeSum<Node>(parent) - safeSum<Node>(node)) + indexOfImpl(parent);
  }

  public:
  /** Gets index of an item in the tree */
  size_t indexOf(const Iterator iter) const noexcept {
    SharedConstNode node = iter.ptr;
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

        n1 -> weight = n3 -> weight ==  1 ? -1 : 0;
        n2 -> weight = n3 -> weight == -1 ?  1 : 0;
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

        n1 -> weight = n3 -> weight == -1 ?  1 : 0;
        n2 -> weight = n3 -> weight ==  1 ? -1 : 0;
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
      return node;
    }
  }



  private:
  /** Inserts a new item into the tree, key mas be unique */
  pair<bool, SharedNode> insertImpl(const Key& key, SharedNode node) noexcept {

    const bool isLeft = cmp(key, node -> key);
    const SharedNode child = isLeft ? node -> left : node -> right;

    if (child == nullptr) {
      SharedNode toInsert = make_shared<Node>(key);
      setChild(node, toInsert, isLeft);
      COUT << "Incing   " << *node << endl;
      node -> weight += isLeft ? -1 : 1;
      node -> reSum();
      return {node -> weight != 0, toInsert};
    } else {
      const auto [incLevel, inserted] = insertImpl(key, child);
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
      return {node -> weight != 0 && incLevel, inserted};
    }
  }

  public:
  /** Inserts a new item into the tree, key mas be unique */
  Iterator insert(const Key& key) noexcept {
    COUT << "Insert   " << key << endl;
    const SharedNode n = root;

    if (n == nullptr) {
      // insert root
      root = make_shared<Node>(key);
      return Iterator{root};
    } else {
      // insert new node
      return Iterator{insertImpl(key, n).second};
    }
  }



  private:
  /** Propagates info about node deletion to it's parents */
  void sendDecreaseUp(SharedNode node, const bool fromLeft, const bool sumOnly) noexcept {
    if (node == nullptr) return;

    const SharedNode parent = node -> parent.lock();
    node -> reSum();

    if (sumOnly) {
      COUT << "Updated  " << *node << endl;
      return sendDecreaseUp(parent, false, true);
    }

    // add node decrease info, inverse of insert
    node -> weight += -1 * (fromLeft ? -1 : 1);

    COUT << "Decrease " << *node << endl;
    if (abs(node -> weight) == 2) {
      COUT << "Rotating " << *node << endl;
      node = rotateNode(node);
      COUT << "Rot res  " << *node << endl;

      if (node -> weight == 0)
        sendDecreaseUp(parent, node -> isLeftChild(), false);
      else
        sendDecreaseUp(parent, false, true);
    } 
    else if (node -> weight == 0)
      sendDecreaseUp(parent, node -> isLeftChild(), false);
    else
      // no further weights updates, just resuming
      sendDecreaseUp(parent, false, true);
  }

  void switchNodes(const SharedNode s, const SharedNode d) noexcept {

    // if d is parent of s -> switch order and handle below
    if (s -> parent.lock() == d) return switchNodes(d, s);

    const SharedNode sl = s -> left;
    const SharedNode sr = s -> right;
    const SharedNode sp = s -> parent.lock();
    const SharedNode dl = d -> left;
    const SharedNode dr = d -> right;
    const SharedNode dp = d -> parent.lock();
    const bool isSLeft = s -> isLeftChild();
    const bool isDLeft = d -> isLeftChild();

    if (s == d -> parent.lock()) {
      // if s is parent of d
      setChild(sp, d, isSLeft);

      if (!isDLeft)
        setChild(d, sl, true);
      else
        setChild(d, sr, false);

      setChild(d, s, isDLeft);
      setChild(s, dl, true);
      setChild(s, dr, false);
    } else {
      // if the nodes have no common edge
      setChild(sp, d, isSLeft);
      setChild(d, sl, true);
      setChild(d, sr, false);

      setChild(dp, s, isDLeft);
      setChild(s, dl, true);
      setChild(s, dr, false);
    }
    swap(s -> weight, d -> weight);
    sendDecreaseUp(s, false, true); // just resum
    sendDecreaseUp(d, false, true); // just resum
  }

  public:
  // 
  void deleteNode(Iterator& iter) noexcept {
    const SharedNode node = iter.unsafePtr();
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
      Iterator sucIter = iter;
      ++sucIter;
      COUT << "Switch   " << sucIter -> key << endl;
      switchNodes(iter.unsafePtr(), sucIter.unsafePtr());
      COUT << "S res it " << *(iter.ptr) << endl;
      COUT << "S res sc " << *(sucIter.ptr) << endl;
      deleteNode(iter);
    } else {
      // Delete node with 1 child
      const SharedNode child = node -> left != nullptr ? node -> left : node -> right;

      setChild(parent, child, isLeft);
      sendDecreaseUp(parent, isLeft, false);
    }
  }
  


  private:
  /** Finds the lower bound for the key given 
   * Does not support equality! */
  SharedConstNode lowerBoundImpl(const Key& key, const SharedConstNode node) const noexcept {
     if (node == nullptr) return nullptr;
     if (cmp(key, node -> key)) {
        COUT << "Lower L  " << *node << endl;
        const auto res = lowerBoundImpl(key, node -> left);
        return res == nullptr ? node : res;
     } else {
        COUT << "Lower R  " << *node << endl;
        return lowerBoundImpl(key, node -> right);
     }
  }

  public:
  /** Finds the lower bound for the key given */
  Iterator lowerBound(const Key& key) const noexcept {
    COUT << "Lower    " << key << endl;
    const Key fake = Key(key.items, -1, key.data);
    return Iterator{lowerBoundImpl(fake, root)};
  }
  /** Finds the lower bound for the key given */
  Iterator lowerBound(const Iterator itr) const noexcept {
    return lowerBound(itr -> key);
  }



  private:
  /** Finds the item before upper bound */
  // TODO reimplement to support correct upper bound
  SharedConstNode upperBoundImpl(const Key& key, const SharedConstNode node) const noexcept {
     if (node == nullptr) return nullptr;
     if (cmp(node -> key, key)) {
        COUT << "Upper R  " << *node << endl;
        const auto res = upperBoundImpl(key, node -> right);
        return res == nullptr ? node : res;
     } else {
        COUT << "Upper L  " << *node << endl;
        return upperBoundImpl(key, node -> left);
     }
  }

  public:
  /** Finds the item before upper bound */
  Iterator upperBound(const Key& key) const noexcept {
    COUT << "Upper    " << key << endl;
    const Key fake = Key(key.items, 1, key.data);
    return Iterator{upperBoundImpl(fake, root)};
  }
  /** Finds the item before upper bound */
  Iterator upperBound(const Iterator itr) const noexcept {
    return upperBound(itr -> key);
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

  private:
  size_t sumBranch(const Iterator low, const Iterator top, const bool fromLeft) const noexcept {
    size_t sum = -1 * safeProdSum<Node>(fromLeft ? low -> left : low -> right)
      + safeProdSum<Node>(low.ptr);

    SharedConstNode node = low.ptr;

    while(true) {
      if (node == top.ptr)
        return sum - safeProdSum<Node>(fromLeft ? node -> right : node -> left);

      const SharedConstNode parent = node -> parent.lock();
      const bool isLeftChild = node -> isLeftChild();
      if (isLeftChild == fromLeft)
        sum += safeProdSum<Node>(parent) - safeProdSum<Node>(node);

      node = parent;
    }
  }

  public:
  size_t sumRange(const size_t i0, const size_t i1) const noexcept {
    const Iterator from = byIndex(i0);
    const Iterator to   = byIndex(i1);
    const Iterator common = commonAncestor(from, to);

    return sumBranch (from, common, true) 
      + sumBranch (to, common, false) 
      - common -> key.items;
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

  private:
  void checkParentChildPointers() const {
    bool rootFound = false;
    size_t nodesFound = 0;
    for (auto itr = begin(); itr != end(); ++itr) {
      const SharedConstNode node = itr.ptr;
      const SharedConstNode parent = node -> parent.lock();
      const SharedConstNode left  = node -> left;
      const SharedConstNode right = node -> right;

      if (parent == nullptr) {
        if (rootFound)
          throw runtime_error("Second root found!");
        rootFound = true;
      } else {
        if (!(parent -> left == node || parent -> right == node)) {
          COUT << "Child/Parent mismatch for " << node -> key.data << endl;
          throw runtime_error("Parent/Child mismatch");
        }
      }
      if (left != nullptr && left -> parent.lock() != node) {
          COUT << "Parent/Child mismatch for " << node -> key.data << endl;
          throw runtime_error("Parent/Child mismatch");
      }
      if (right != nullptr && right -> parent.lock() != node) {
          COUT << "Parent/Child mismatch for " << node -> key.data << endl;
          throw runtime_error("Parent/Child mismatch");
      }

      ++nodesFound;
    }
    if (nodesFound != size()) {
          COUT << "Number of node/Sum mismatch " << nodesFound << " " << size() << endl;
          throw runtime_error("Number of node/Sum mismatch");
    }
  }

  size_t weightTest(const SharedConstNode node) const {
    if (node == nullptr) return 0;
    const ssize_t l = weightTest(node -> left);
    const ssize_t r = weightTest(node -> right);
    const ssize_t sum =  r - l;
    if (sum != node -> weight) {
      COUT << "Weight assertion failed for " << node -> key << ", has " << (int)node -> weight << ", got " << sum << endl;
      throw runtime_error("Weight assertion failed");
    }
    return 1 + (l > r ? l : r);
  }

  size_t sumsTest(const SharedConstNode node) const {
    if (node == nullptr) return 0;
    const ssize_t l = sumsTest(node -> left);
    const ssize_t r = sumsTest(node -> right);
    const size_t sum =  l + r + 1;
    if (sum != node -> sum) {
      COUT << "Sum assertion failed for " << node -> key << endl;
      throw runtime_error("Sum assertion failed");
    }
    return sum;
  }

  size_t prodSumsTest(const SharedConstNode node) const {
    if (node == nullptr) return 0;
    const ssize_t l = prodSumsTest(node -> left);
    const ssize_t r = prodSumsTest(node -> right);
    const size_t sum =  l + r + node -> key.items;
    if (sum != node -> sellSum) {
      COUT << "Sell sum assertion failed for " << node -> key << endl;
      throw runtime_error("Sell sum assertion failed");
    }
    return sum;
  }

  void compareTest(const SharedConstNode node) const {
      if (node == nullptr) return;
      const SharedConstNode left  = node -> left;
      const SharedConstNode right = node -> right;

      if (left != nullptr) {
        if(!cmp(left -> key, node -> key)) {
          COUT << "Compare assertion failed for " << node -> key << endl;
          throw runtime_error("Compare assertion failed");
        }
        compareTest(left);
      }

      if (right != nullptr) {
        if(!cmp(node -> key, right -> key)) {
          COUT << "Compare assertion failed for " << node -> key << endl;
          throw runtime_error("Compare assertion failed");
        }
        compareTest(right);
      }
  }

  public:
  void validate() const {
    checkParentChildPointers();
    sumsTest(root);
    prodSumsTest(root);
    compareTest(root);
    weightTest(root);
  }

  template<typename Product>
  friend struct Bestsellers;
};

template<typename Product>
struct Bestsellers {
  private:
  using PHolder = Holder<Product>;
  using Tree = AVLTree<PHolder, BaseNode<Product>, HolderCmp<PHolder>>;
  using Map = unordered_map<Product, typename Tree::Iterator>;

  Tree tree;
  Map map;

  void throwIndexOutOf(size_t index) const {
    throw std::out_of_range("Index is out of range: "s + to_string(index));
  }

  public:
  // The total number of tracked products
  size_t products() const { return map.size(); }

  void sell(const Product& p, size_t amount) {
    const auto res = map.find(p);
    if (res == map.end()) {
      const auto itr = tree.insert(PHolder(amount, p));
      map.emplace(p, itr);
    } else {
      auto prev = res -> second;
      const size_t prevSum = prev -> key.items;
      tree.deleteNode(prev);
      const auto itr = tree.insert(PHolder(amount + prevSum, p));
      map.insert_or_assign(p, itr);
    }
  }

  // The most sold product has rank 1
  size_t rank(const Product& p) const {
    const auto res = map.find(p);
    if (res == map.end())
      throw std::out_of_range("Item not found");

    return 1 + tree.indexOf(res -> second);
  }
  const Product& product(size_t rank) const {
    if (--rank >= map.size()) throwIndexOutOf(rank);
    return tree.byIndex(rank) -> key.data;
  }

  // How many copies of product with given rank were sold
  size_t sold(size_t rank) const { return sold(rank, rank); }
  // The same but sum over interval of products (including from and to)
  // It must hold: sold(x) == sold(x, x)
  size_t sold(size_t from, size_t to) const {
    if (--from >= map.size()) throwIndexOutOf(from);
    if (--to   >= map.size()) throwIndexOutOf(to);
    if (from > to        )
      throw std::out_of_range("From is greater than to: "s + to_string(from) + " > " + to_string(to));

    return tree.sumRange(from, to);
  }

  // The smallest (resp. largest) rank with sold(rank) == sold(r)
  size_t first_same(size_t rank) const {
    if (--rank >= map.size()) throwIndexOutOf(rank);

    return 1 + tree.indexOf(tree.lowerBound(tree.byIndex(rank)));
  }

  size_t last_same(size_t rank) const {
    if (--rank >= map.size()) throwIndexOutOf(rank);

    return 1 + tree.indexOf(tree.upperBound(tree.byIndex(rank)));
  }

  void print() const noexcept { tree.print(); }
  void validate() const { tree.validate(); }
};

#ifndef __PROGTEST__

void test1() {
  Bestsellers<std::string> T;
   T.sell("coke", 32);
   T.sell("bread", 1);
   assert(T.products() == 2);
   T.sell("ham", 2);
   T.sell("mushrooms", 12);
 
   assert(T.products() == 4);
   assert(T.rank("ham") == 3);
   assert(T.rank("coke") == 1);
   assert(T.sold(1, 3) == 46);
   assert(T.product(2) == "mushrooms");
 
   T.sell("ham", 11);
   assert(T.products() == 4);
   assert(T.product(2) == "ham");
   assert(T.sold(2) == 13);
   assert(T.sold(2, 2) == 13);
   assert(T.sold(1, 2) == 45);
}
 
void test2() {
# define CATCH(expr) \
  try { expr; assert(0); } catch (const std::out_of_range&) { assert(1); };

  Bestsellers<std::string> T;
  T.sell("coke", 32);
  T.sell("bread", 1);

  CATCH(T.rank("ham"));
  CATCH(T.product(3));
  CATCH(T.sold(0));
  CATCH(T.sold(9));
  CATCH(T.sold(0, 1));
  CATCH(T.sold(3, 2));
  CATCH(T.sold(1, 9));

#undef CATCH
}

struct Muhaha {
  bool operator<(const Muhaha& other) const noexcept { return true; }
  bool operator==(const Muhaha& other) const noexcept { return true; }
};

namespace std {
  template <> struct hash<Muhaha> {
    size_t operator()(const Muhaha & x) const { return 0; }
  };
}

void testRankRanges() {
  Bestsellers<ssize_t> T;
  T.sell(0, 10);
  T.sell(1, 10);
  T.sell(2, 10);
  T.sell(3, 20);
  T.sell(4, 20);
  T.sell(5, 20);
  T.sell(6, 30);
  T.sell(7, 30);
  T.sell(8, 30);
  assert(T.first_same(1) == 1);
  assert(T.first_same(2) == 1);
  assert(T.first_same(3) == 1);
  assert(T.last_same(1) == 3);
  assert(T.last_same(2) == 3);
  assert(T.last_same(3) == 3);
  assert(T.first_same(4) == 4);
  assert(T.first_same(5) == 4);
  assert(T.first_same(6) == 4);
  assert(T.last_same(4) == 6);
  assert(T.last_same(5) == 6);
  assert(T.last_same(6) == 6);
  assert(T.first_same(7) == 7);
  assert(T.first_same(8) == 7);
  assert(T.first_same(9) == 7);
  assert(T.last_same(7) == 9);
  assert(T.last_same(8) == 9);
  assert(T.last_same(9) == 9);
}

void testCompile() {
  // Bestsellers<Muhaha> b;
  // b.sell(Muhaha(), 1);
  // b.rank(Muhaha());
  // b.product(1);
  // b.sold(1);
  // b.sold(1, 1);
  // b.first_same(1);
  // b.last_same(1);
}
size_t random(const size_t min, const size_t max) {
  return rand() % (max - min + 1) + min;
}
size_t random(const size_t max) { return random(0, max); }

void bulkTest() {
  if constexpr (true) srand (time(NULL));

  for (size_t t = 0; t < 4; ++t) {
    Bestsellers<ssize_t> b;
    const size_t maxName = 1 << 10;
    const size_t maxCount = 1 << 10;
    for (size_t i = 0; i < 1 << 16; ++i) {
      b.sell(random(maxName), random(maxCount));
      // b.print();
      b.validate();

      const size_t rank = random(1, b.products() / 2 + 1);
      assert(rank == b.rank(b.product(rank)));
      const size_t rank2 = random(b.products() / 2 + 1, b.products());
      if (rank < rank2)
        b.sold(rank, rank2);
      else b.sold(rank2, rank);

      b.first_same(rank);
      b.last_same(rank);
    }
  }
}

int main() {
  bulkTest();
  test1();
  test2();
  testRankRanges();
  // testCompile();
  return 0;
}

#endif
