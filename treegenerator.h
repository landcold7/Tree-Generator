/*

MIT License

Copyright (c) 2019 ouuan

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

// See https://github.com/ouuan/Tree-Generator to get latest version or bug
// tracker.

#ifndef TREE_GENERATOR_BY_OUUAN_
#define TREE_GENERATOR_BY_OUUAN_ 1

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdlib>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
#include <iostream>
#include <queue>
#include <random>
#include <string>
#include <vector>

#include "../../bits/debug.h"

namespace tree_generator_by_ouuan {

typedef std::pair<int, int> pii;

std::mt19937 Next(std::chrono::steady_clock::now().time_since_epoch().count());

int DefaultRandmInt(int l, int r) {
  int out = Next() % (r - l + 1) + l;
  return out >= l ? out : out + r - l + 1;
}

int (*randint)(int, int) = DefaultRandmInt;

void DefaultOutPutEdge(std::ostream& os, int u, int parent) {
  if (randint(0, 1) == 0) {
    os << u + 1 << ' ' << parent + 1 << std::endl;
  }
  else {
    os << parent + 1 << ' ' << u + 1 << std::endl;
  }
}

void (*OutputEdge)(std::ostream&, int, int) = DefaultOutPutEdge;

class Tree {
 private:
  std::vector<int> p, id, eid;

 public:
  Tree() {
    p.push_back(-1);
    id.push_back(0);
  }

  Tree(int n) {
    assert(n > 0);
    p.push_back(-1);
    id.push_back(0);
    if (n > 1) {
      Random(n - 1, 0);
    }
  }

  Tree(const std::string& s) {
    p.push_back(-1);
    id.push_back(0);

    std::function<size_t(const std::string&, size_t)> FindComma =
        [](const std::string& str, size_t ptr) {
          while (ptr < str.size() && str[ptr] != ',') ++ptr;
          return ptr;
        };

    std::function<size_t(const std::string&, size_t)> FindLetter =
        [](const std::string& str, size_t ptr) {
          while (ptr < str.size() && !isalpha(str[ptr])) ++ptr;
          return ptr;
        };

    std::function<void(const std::string&)> Error =
        [](const std::string& func) {
          std::cerr << "Error: the number of parentrameters for " << func
                    << " is not correct." << std::endl;
          exit(1);
        };

    size_t pos = 0;

    while (pos < s.size()) {
      if (pos + 1 >= s.size()) {
        std::cerr << "Error: can't recognize the tree type abbreviation, it's "
                     "too short.\n";
        exit(1);
      }
      std::string type = s.substr(pos, 2);
      pos += 2;
      for_each(type.begin(), type.end(), [](char& x) { x = tolower(x); });
      int nextLetter = FindLetter(s, pos);
      if (type == "lh") {
        int nextComma = FindComma(s, pos);
        assert(nextComma < nextLetter);
        int n = atoi(s.substr(pos, nextComma - pos).c_str());
        pos = nextComma + 1;
        nextComma = FindComma(s, pos);
        assert(nextComma < nextLetter);
        double low, high;
        sscanf(s.substr(pos, nextComma - pos).c_str(), "%lf", &low);
        pos = nextComma + 1;
        nextComma = FindComma(s, pos);
        assert(nextComma < nextLetter);
        sscanf(s.substr(pos, nextComma - pos).c_str(), "%lf", &high);
        pos = nextComma + 1;
        nextComma = FindComma(s, pos);
        assert(nextComma >= nextLetter);
        int parent = atoi(s.substr(pos, nextComma - pos).c_str());
        pos = nextLetter;
        LowHigh(n, low, high, parent);
        continue;
      }

      if (type == "lm") {
        int nextComma = FindComma(s, pos);
        assert(nextComma < nextLetter);
        int n = atoi(s.substr(pos, nextComma - pos).c_str());
        pos = nextComma + 1;
        nextComma = FindComma(s, pos);
        assert(nextComma < nextLetter);
        int k = atoi(s.substr(pos, nextComma - pos).c_str());
        pos = nextComma + 1;
        nextComma = FindComma(s, pos);
        assert(nextComma < nextLetter);
        double low, high;
        sscanf(s.substr(pos, nextComma - pos).c_str(), "%lf", &low);
        pos = nextComma + 1;
        nextComma = FindComma(s, pos);
        assert(nextComma < nextLetter);
        sscanf(s.substr(pos, nextComma - pos).c_str(), "%lf", &high);
        pos = nextComma + 1;
        nextComma = FindComma(s, pos);
        assert(nextComma >= nextLetter);
        int parent = atoi(s.substr(pos, nextComma - pos).c_str());
        pos = nextLetter;
        LowHighMaxDegree(n, k, low, high, parent);
        continue;
      }

      std::vector<int> parentr;
      while (1) {
        int nextComma = FindComma(s, pos);
        parentr.push_back(atoi(s.substr(pos, nextComma - pos).c_str()));
        pos = nextComma + 1;
        if (nextComma >= nextLetter) {
          pos = nextLetter;
          break;
        }
      }
      if (type == "nd") {
        if (parentr.size() == 1)
          AddNode(parentr[0]);
        else
          Error("AddNode");
      }
      else if (type == "rd") {
        if (parentr.size() == 2)
          Random(parentr[0], parentr[1]);
        else
          Error("Random");
      }
      else if (type == "tl") {
        if (parentr.size() == 3)
          Tall(parentr[0], parentr[1], parentr[2]);
        else
          Error("Tall");
      }
      else if (type == "ch") {
        if (parentr.size() == 2)
          Chain(parentr[0], parentr[1]);
        else
          Error("Chain");
      }
      else if (type == "st") {
        if (parentr.size() == 2)
          Star(parentr[0], parentr[1]);
        else
          Error("Star");
      }
      else if (type == "fl") {
        if (parentr.size() == 2)
          Flower(parentr[0], parentr[1]);
        else
          Error("Flower");
      }
      else if (type == "md") {
        if (parentr.size() == 3)
          MaxDegree(parentr[0], parentr[1], parentr[2]);
        else
          Error("MaxDegree");
      }
      else if (type == "cp") {
        if (parentr.size() == 3)
          Complete(parentr[0], parentr[1], parentr[2]);
        else
          Error("Complete");
      }
      else if (type == "bi") {
        if (parentr.size() == 2)
          Binary(parentr[0], parentr[1]);
        else
          Error("Binary");
      }
      else if (type == "cb") {
        if (parentr.size() == 2)
          CompleteBinary(parentr[0], parentr[1]);
        else
          Error("CompleteBinary");
      }
      else if (type == "sw") {
        if (parentr.size() == 2)
          Silkworm(parentr[0], parentr[1]);
        else
          Error("Silkworm");
      }
      else if (type == "al") {
        if (parentr.size() == 3)
          AddLeaves(parentr[0], parentr[1], parentr[2]);
        else
          Error("AddLeaves");
      }
      else {
        std::cerr << "Error: can't recognize the tree type abbreviation "
                  << type << "." << std::endl;
        exit(1);
      }
    }
  }

  int Size() const { return id.size(); }

  void AddNode(int par) {
    assert(par >= 0);
    assert(par < Size());
    // New node that will be added.
    int u = id.size();
    // u <---> parent
    id.push_back(u);
    p.push_back(par);
    eid.push_back(u);
    // trace(id, p, eid);
  }

  void Random(int n, int parent) {
    int sz = Size();
    assert(n > 0);
    assert(parent >= 0);
    assert(parent < sz);
    AddNode(parent);
    if (n == 1) return;
    if (n == 2) {
      AddNode(sz);
      return;
    }
    // https://en.wikipedia.org/wiki/Pr%C3%BCfer_sequence
    std::vector<int> prufer, cnt;
    std::vector<std::vector<int> > g;
    g.resize(n);
    cnt.resize(n, 0);
    for (int i = 0; i < n - 2; ++i) {
      int x = randint(0, n - 1);
      prufer.push_back(x);
      ++cnt[x];
    }
    std::priority_queue<int> q;
    for (int i = 0; i < n; ++i) {
      if (!cnt[i]) {
        q.push(i);
      }
    }
    trace(n, prufer, cnt, q);
    for (auto v : prufer) {
      int u = q.top();
      g[u].push_back(v);
      g[v].push_back(u);
      // trace(v, u);
      q.pop();
      if (--cnt[v] == 0) q.push(v);
    }
    int x = q.top();
    q.pop();
    int y = q.top();
    g[x].push_back(y);
    g[y].push_back(x);
    // trace(x, y);

    std::queue<int> bfs;
    bfs.push(0);
    int root = sz;
    while (!bfs.empty()) {
      int u = bfs.front();
      cnt[u] = 1;
      bfs.pop();
      for (auto v : g[u]) {
        if (cnt[v] == 0) {
          AddNode(root);
          bfs.push(v);
        }
      }
      ++root;
    }
  }

  void LowHigh(int n, double low, double high, int parent) {
    int sz = Size();
    assert(n > 0);
    assert(low >= 0);
    assert(high <= 1);
    assert(high >= low);
    assert(parent >= 0);
    assert(parent < sz);
    AddNode(parent);
    for (int i = 1; i < n; ++i) {
      AddNode(randint(round((i - 1) * low), round((i - 1) * high)) + sz);
    }
  }

  void Tall(int n, int k, int parent) {
    int sz = Size();
    assert(n > 0);
    assert(k > 0);
    assert(parent >= 0);
    assert(parent < sz);
    AddNode(parent);
    for (int i = sz + 1; i < sz + n; ++i) {
      AddNode(randint(std::max(sz, i - k), i - 1));
    }
  }

  void Chain(int n, int parent) {
    assert(n > 0);
    assert(parent >= 0);
    assert(parent < Size());
    Tall(n, 1, parent);
  }

  void Star(int n, int parent) {
    int sz = Size();
    assert(n > 0);
    assert(parent >= 0);
    assert(parent < sz);
    AddNode(parent);
    for (int i = sz + 1; i < sz + n; ++i) {
      AddNode(sz);
    }
  }

  void Flower(int n, int parent) {
    assert(n > 0);
    assert(parent >= 0);
    assert(parent < Size());
    Star(n, parent);
  }

  void MaxDegree(int n, int k, int parent) {
    int sz = Size();
    assert(n > 0);
    assert(k >= 2);
    assert(parent >= 0);
    assert(parent < sz);
    AddNode(parent);
    __gnu_pbds::tree<pii, __gnu_pbds::null_type, std::less<pii>,
                     __gnu_pbds::rb_tree_tag,
                     __gnu_pbds::tree_order_statistics_node_update>
        remain;
    remain.insert(pii(sz, k - 1));
    for (int i = sz + 1; i < sz + n; ++i) {
      auto it = remain.find_by_order(randint(0, remain.size() - 1));
      int u = it->first;
      int t = it->second;
      remain.erase(it);
      if (t > 1) remain.insert(pii(u, t - 1));
      AddNode(u);
      remain.insert(pii(i, k - 1));
    }
  }

  void LowHighMaxDegree(int n, int k, double low, double high, int parent) {
    int sz = Size();
    assert(n > 0);
    assert(k >= 2);
    assert(low >= 0);
    assert(high <= 1);
    assert(parent >= 0);
    assert(parent < sz);
    AddNode(parent);
    __gnu_pbds::tree<pii, __gnu_pbds::null_type, std::less<pii>,
                     __gnu_pbds::rb_tree_tag,
                     __gnu_pbds::tree_order_statistics_node_update>
        remain;
    remain.insert(pii(sz, k - 1));
    for (int i = sz + 1; i < sz + n; ++i) {
      auto it = remain.find_by_order(randint(
          round((remain.size() - 1) * low), round((remain.size() - 1) * high)));
      int u = it->first;
      int t = it->second;
      remain.erase(it);
      if (t > 1) remain.insert(pii(u, t - 1));
      AddNode(u);
      remain.insert(pii(i, k - 1));
    }
  }

  void Complete(int n, int k, int parent) {
    int sz = Size();
    assert(n > 0);
    assert(k >= 2);
    assert(parent >= 0);
    assert(parent < sz);
    AddNode(parent);
    for (int i = sz + 1; i < sz + n; ++i)
      AddNode(sz + ceil(1.0 * (i - sz) / (k - 1) - 1e-9) - 1);
  }

  void Binary(int n, int parent) {
    assert(n > 0);
    assert(parent >= 0);
    assert(parent < Size());
    MaxDegree(n, 3, parent);
  }

  void CompleteBinary(int n, int parent) {
    assert(n > 0);
    assert(parent >= 0);
    assert(parent < Size());
    Complete(n, 3, parent);
  }

  void Silkworm(int n, int parent) {
    int sz = Size();
    assert(n > 0);
    assert(parent >= 0);
    assert(parent < sz);
    int chain_len = (n + 1) / 2;
    Chain(chain_len, parent);
    for (int i = sz; i + chain_len < sz + n; ++i) {
      AddNode(i);
    }
  }

  void AddLeaves(int n, int l, int r) {
    assert(n > 0);
    assert(l >= 0);
    assert(r < Size());
    assert(l <= r);
    for (int i = 1; i <= n; ++i) {
      AddNode(randint(l, r));
    }
  }

  void shuffleNodes(int from = 0) {
    for (int i = 0; i < from; ++i) id[i] = i;
    for (size_t i = from; i < id.size(); ++i) {
      id[i] = i;
      std::swap(id[i], id[randint(from, i)]);
    }
  }

  void shuffleEdges() {
    for (size_t i = 0; i < eid.size(); ++i) {
      eid[i] = i + 1;
      std::swap(eid[i], eid[randint(0, i)]);
    }
  }

  void Resize(int n) {
    assert(n > 0);
    int sz = Size();
    if (sz < n)
      AddLeaves(n - sz, 0, sz - 1);
    else if (sz > n) {
      p.resize(n);
      id.resize(n);
      eid.resize(n - 1);
      for (int i = 0; i < n; ++i) id[i] = i;
      for (int i = 0; i < n - 1; ++i) eid[i] = i + 1;
    }
  }

  void PrintEdge(int x, std::ostream& os = std::cout) const {
    OutputEdge(os, id[eid[x]], id[p[eid[x]]]);
  }

  int Parent(int u) const { return p[u]; }

  friend std::ostream& operator<<(std::ostream& os, const Tree& t) {
    for (int i = 0; i < (int)t.eid.size(); ++i) {
      t.PrintEdge(i, os);
    }
    return os;
  }
};
}  // namespace tree_generator_by_ouuan

#endif
