#ifndef __PROGTEST__
#include <cassert>
#include <cstdint>
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
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <stack>
#include <queue>
#include <random>

using ChristmasTree = size_t;

struct TreeProblem {
  int max_group_size;
  std::vector<uint64_t> gifts;
  std::vector<std::pair<ChristmasTree, ChristmasTree>> connections;
};

#endif

using namespace std;

using Args = TreeProblem;
using Point = ChristmasTree;
using Map = vector<vector<Point>>;
using BackMap = vector<Point>;
using StartingList = vector<Point>;
using Cache = vector<size_t>;

struct Tree {
  const Map map;
  const BackMap back;
  const StartingList starting;
};


const Point start = 0;

template<typename T>
ostream& operator<<(ostream& out, const vector<T>& data) {
  bool isFirst = true;
  out << "(";
  for (const auto& target: data) {
    if (isFirst) isFirst = false;
    else out << ",";

    out << target;
  }
  return out << ")";
}

void printCache(const Cache& cache, size_t rowWidth, ostream& out = cout) {
  for (size_t i = 0; i < cache.size() / rowWidth; ++i) {
    bool isFirst = true;
    for (size_t j = 0; j < rowWidth; ++j) {
      if (isFirst) isFirst = false;
      else out << ",";

      size_t t = cache[i * rowWidth + j];
      if (t == 0) out << "*"; else out << t - 1;
    }
    out << "\n";
  }
}

inline Map buildVerticies(const Args& args) {
  Map map(args.gifts.size());

  for (const auto& conn : args.connections) {
    map[conn.first].emplace_back(conn.second);
    map[conn.second].emplace_back(conn.first);
  }

  return map;
}

Tree buildTree(const Args& args, const Map& vert) {
  const size_t totalItems = args.gifts.size();
  Map map(totalItems);
  BackMap back(totalItems);
  StartingList starting;

  queue<pair<Point, Point>> q;
  const size_t invalidState = (Point) -1;
  q.emplace(make_pair(start, invalidState));

  while (!q.empty()) {
    const auto [item, parent] = q.front();
    q.pop();

    for (const auto curr : vert[item]) {
      if (curr == parent) continue;

      map[item].emplace_back(curr);
      back[curr] = item;
      q.emplace(make_pair(curr, item));
    }

    if (map[item].empty()) {
      starting.emplace_back(item);
    }
  }

  return Tree{ move(map), move(back), move(starting) };
}

size_t resolveNodeIterativeTopSort(
    const Args& args,
    const Tree& tree
) {
  const size_t totalItems = args.gifts.size();
  const bool twoAllowed = args.max_group_size == 2;

  const Map& map = tree.map;
  const BackMap& back = tree.back;
  const StartingList& starting = tree.starting;
  const vector<uint64_t>& gifts = args.gifts;

  Cache cache((1 + args.max_group_size) * totalItems);

  queue<Point> q;

  vector<uint32_t> topSortArray(totalItems); 

  for (const auto& edge : back) {
    ++topSortArray[edge];
  }
  --topSortArray[0];
  // cout << "Top " << topSortArray << endl;

  for (const Point p : starting) {
    q.emplace(p);
  }

  while (!q.empty()) {
    const Point item = q.front();
    q.pop();

    // cout << "Handling " << item << endl;

    size_t sumNone = 0; // not guarded
    size_t sumOne  = gifts[item]; // guarded by one or two
    ssize_t maxDiff = 0;

    // compute value
    for (const auto& curr : map[item]) {

      const ssize_t byNone = cache[curr] - 1;
      const ssize_t byOne  = cache[curr + totalItems] - 1;

      if (twoAllowed) {
        const ssize_t byTwo  = cache[curr + 2 * totalItems] - 1;
        sumNone += max(max(byNone, byOne), byTwo);
        sumOne  += byNone;
        maxDiff = max(maxDiff, byTwo - byNone);
      } else {
        sumNone += max(byNone, byOne);
        sumOne  += byNone;
      }
    }

    cache[item] = sumNone + 1;
    cache[item + totalItems] = sumOne + maxDiff + 1;
    if (twoAllowed) {
      cache[item + 2 * totalItems] = sumOne + 1;
    }

    {
      // schedule next ones
      const Point curr = back[item];
      --topSortArray[curr];
      if (topSortArray[curr] == 0) {
        q.emplace(curr);
      }
    }

    // printCache(cache, totalItems);
  }

  if (twoAllowed) {
    return max(
        max(cache[start], cache[start + totalItems]),
        cache[start + 2 * totalItems]
        ) - 1;
  } else {
    return max(cache[start], cache[start + totalItems]) - 1;
  }
}

uint64_t solve(const Args& args) {
  if (args.gifts.size() == 0)
    return 0;

  const Tree tree = buildTree(args, buildVerticies(args));
  // cout << tree.map << endl;
  // cout << tree.back << endl;
  // cout << tree.starting << endl;
  return resolveNodeIterativeTopSort(args, tree);
}

#ifndef __PROGTEST__

using TestCase = std::pair<uint64_t, TreeProblem>;

const std::vector<TestCase> BASIC_TESTS = {
  { 3, { 1, { 1, 1, 1, 2 }, { {0,3}, {1,3}, {2,3} }}},
  { 4, { 1, { 1, 1, 1, 4 }, { {0,3}, {1,3}, {2,3} }}},
  { 57, { 1, {
    17, 11, 5, 13, 8, 12, 7, 4, 2, 8,
  }, {
    {1, 4}, {6, 1}, {2, 1}, {3, 8}, {8, 0}, {6, 0}, {5, 6}, {7, 2}, {0, 9},
  }}},
  { 85, { 1, {
    10, 16, 13, 4, 19, 8, 18, 17, 18, 19, 10,
  }, {
    {9, 7}, {9, 6}, {10, 4}, {4, 9}, {7, 1}, {0, 2}, {9, 2}, {3, 8}, {2, 3}, {5, 4},
  }}},
  { 79, { 1, {
    8, 14, 11, 8, 1, 13, 9, 14, 15, 12, 1, 11,
  }, {
    {9, 1}, {1, 2}, {1, 4}, {5, 10}, {7, 8}, {3, 7}, {11, 3}, {11, 10}, {6, 8}, {0, 1}, {0, 3},
  }}},
  { 102, { 1, {
    15, 10, 18, 18, 3, 4, 18, 12, 6, 19, 9, 19, 10,
  }, {
    {10, 2}, {11, 10}, {6, 3}, {10, 8}, {5, 3}, {11, 1}, {9, 5}, {0, 4}, {12, 3}, {9, 7}, {11, 9}, {4, 12},
  }}},
  { 93, { 1, {
    1, 7, 6, 18, 15, 2, 14, 15, 18, 8, 15, 1, 5, 6,
  }, {
    {0, 13}, {6, 12}, {0, 12}, {7, 8}, {8, 3}, {12, 11}, {12, 1}, {10, 12}, {2, 6}, {6, 9}, {12, 7}, {0, 4}, {0, 5},
  }}},
};

const std::vector<TestCase> BONUS_TESTS = {
  { 3, { 2, { 1, 1, 1, 2 }, { {0,3}, {1,3}, {2,3} }}},
  { 5, { 2, { 1, 1, 1, 4 }, { {0,3}, {1,3}, {2,3} }}},
};

void printLine(ostream& out = cout) {
  for (int i = 0; i < 80; ++i)
    out << "-";
  out << "\n";
}

void test(const std::vector<TestCase>& T) {
  int i = 0;
  for (auto &[s, t] : T) {

    printLine();
    cout << "Checking " << i << endl;
    printLine();

    size_t solution = solve(t);

    if (s != solution) {
      std::cout << "ERROR in " << i << " (mismatch " << s << " != " << solution << ")"<< std::endl;
    } else {
      cout << "PASSED" << endl;
    }

    i++;
  }
  std::cout << "Finished" << std::endl;
}

int main() {
  printLine();
  cout << "BASIC TEST" << endl;
  printLine();
  test(BASIC_TESTS);
  printLine();
  cout << "BONUS TEST" << endl;
  printLine();
  test(BONUS_TESTS);
}

#endif


