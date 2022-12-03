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
using Map = unordered_multimap<Point, Point>;
using Cache = vector<size_t>;


const Point start = 0;

void printMap(const Map& map, ostream& out = cout) {
  bool isFirst = true;
  out << "{";
  for (const auto& vert : map) {

    if (isFirst) isFirst = false;
    else out << ",";

    out << "(" << vert.first << "," << vert.second << ")";
  }
  out << "}" << endl;
}


Map buildVerticies(const Args& args) {
  Map map;
  for (const auto& conn : args.connections) {
    map.emplace(conn);
    map.emplace(make_pair(conn.second, conn.first));
  }
  return map;
}

Map buildTree(const Map& vert) {
  Map out;

  queue<pair<Point, Point>> q;
  q.emplace(make_pair(start, (Point) -1));

  while (!q.empty()) {
    const auto [item, parent] = q.front();
    q.pop();

    const auto range = vert.equal_range(item);
    for (auto itr = range.first; itr != range.second; ++itr) {

      const Point curr = itr -> second;
      if (curr == parent) continue;

      out.emplace(make_pair(item, curr));
      q.emplace(make_pair(curr, item));
    }
  }

  return out;
}

size_t resolveNode(
    const Args& args,
    const Map& tree,
    Cache& cache,
    const Point node,
    const bool shouldGuard
    ) {

  const size_t cacheIndex = node + (shouldGuard ? args.gifts.size() : 0);
  if (cache[cacheIndex] != 0)
    return cache[cacheIndex] - 1;

  size_t max = 0;

  const auto range = tree.equal_range(node);
  for (auto itr = range.first; itr != range.second; ++itr) {
    const Point curr = itr -> second;

    size_t r = resolveNode(args, tree, cache, curr, false);

    if (!shouldGuard) {
      const size_t r2 = resolveNode(args, tree, cache, curr, true);
      if (r2 > r) r = r2;
    }
    max += r;
  }

  const size_t res = max + (shouldGuard ? args.gifts[node] : 0);
  cache[cacheIndex] = res + 1;

  return res;
}

uint64_t solve(const Args& args) {
  if (args.gifts.size() == 0)
    return 0;

  // printMap(buildVerticies(args));
  const Map tree = buildTree(buildVerticies(args));
  Cache cache((1 + args.max_group_size) * args.gifts.size());

  size_t max1 = resolveNode(args, tree, cache, 0, true);
  size_t max2 = resolveNode(args, tree, cache, 0, false);

  return max1 > max2 ? max1 : max2;
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

void test(const std::vector<TestCase>& T) {
  int i = 0;
  for (auto &[s, t] : T) {
    size_t solution = solve(t);
    if (s != solution)
      std::cout << "Error in " << i << " (mismatch " << s << " != " << solution << ")"<< std::endl;
    i++;
  }
  std::cout << "Finished" << std::endl;
}

int main() {
  test(BASIC_TESTS);
  // test(BONUS_TESTS);
}

#endif


