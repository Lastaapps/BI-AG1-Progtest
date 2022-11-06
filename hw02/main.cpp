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

template < typename Key>
class AVLTree {
  struct Node {
    using SharedNode = shared_ptr<Node>;
    const Key key;
    size_t sum = 0;
    SharedNode left = nullptr, right = nullptr;
    weak_ptr<Node> parent;
    char weight = 0;

    Node(const Key& k) : key(k) {}

    bool isLeftChild() const {
      const SharedNode lock = parent.lock();
      // return lock == nullptr || lock.get() == this;
      return lock == nullptr || lock -> left.get() == this;
    }

    friend ostream& operator<<(ostream& out, const Node& node) {
      const SharedNode p = node.parent.lock();
      const string pout = p == nullptr ? "null"s : to_string(p -> key);
      out << "{"
        << node.key
        << ", w: " << (int)node.weight
        << ", s: " << node.sum
        << ", l: " << node.isLeftChild()
        << ", p: " << pout
        << "}";
      return out;
    }
  };
  using SharedNode = shared_ptr<Node>;
  SharedNode root = nullptr;

  struct Iterator {
    SharedNode ptr;
    Node* operator->() {
      return ptr.get();
    }
    Node& operator*() {
      return *ptr;
    }
    void operator++() {
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
    bool operator==(const Iterator& other) {
      return ptr == other.ptr;
    }
    bool operator!=(const Iterator& other) {
      return ptr != other.ptr;
    }
  };

  public:
  Iterator begin() {
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

  Iterator end() {
    return Iterator({nullptr});
  }

  private:

  void setChild(SharedNode parent, SharedNode child, bool isLeft) {
    if (parent != nullptr)
      (isLeft ? parent -> left : parent -> right) = child;
    else
      root = child;

    if (child != nullptr) {
      child -> parent = parent;
    }
  }

  SharedNode rotateNode(SharedNode node) {
    const char weight = node -> weight;
    SharedNode parent = node -> parent.lock();
    const bool isLeft = node -> isLeftChild();

    // TODO sum
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
        return n3;
      }
    } else {
      throw runtime_error("You fucked up!");
    }
  }

  bool insert_node(const Key& key, SharedNode node) {

    bool isLeft = key < node -> key;
    SharedNode child = isLeft ? node -> left : node -> right;

    if (child == nullptr) {
      SharedNode toInsert = make_shared<Node>(key);
      setChild(node, toInsert, isLeft);
      node -> weight += isLeft ? -1 : 1;
      cout << "Incing   " << *node << endl;
      return true;
    } else {
      const bool incLevel = insert_node(key, child);
      node -> sum += 1;
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
  void insert(const Key& key) {
    cout << "Insert   " << key << endl;
    SharedNode n = root;

    if (n == nullptr) {
      root = make_shared<Node>(key);
      return;
    }

    insert_node(key, n);
  }

  void printImpl(ostream& out, SharedNode node, size_t depth) {
    if (node == nullptr) {
      // out << "null" << endl;
      return;
    }
    for (size_t i = 0;  i < depth; ++i) out << " ";
    out << *node << endl;
    printImpl(out, node -> left, depth + 1);
    printImpl(out, node -> right, depth + 1);
  }

  void print(ostream& out = cout) { printImpl(out, root, 0); }
};

template<typename Product>
struct Bestsellers {
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

  // Double right rotation
  {
    AVLTree<ssize_t> tree;
    for (size_t i = 0; i < 20; ++i) {
      tree.insert(i % 2 == 0 ? i / 2 : 19 - i / 2);
      tree.print();
    }
    cout << "\n\n\n";
  }

  // Iterator
  //  AVLTree<ssize_t> tree;
  //  for (size_t i = 20; i > 0; --i) tree.insert(i - 1);
  //  tree.print();
  //   for (auto itr = tree.begin(); itr != tree.end(); ++itr)
  //     cout << *itr << endl;
  //  cout << "\n\n\n";
  //}
}

// void test1() {
//   Bestsellers<std::string> T;
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
// }
// 
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
  // test1();
  // test2();
  return 0;
}

#endif


