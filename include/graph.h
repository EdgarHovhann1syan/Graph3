#ifndef GRAPH_H
#define GRAPH_H
#include <vector>
#include <limits>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <algorithm>
#include <stack>
class Graph 
{
public:
    Graph(std::size_t countOfvertices, bool is_directed = true);
    void add_edge(std::size_t u, std::size_t v, int weight);
    void remove_edge(std::size_t u, std::size_t v);
    void dfs(std::size_t u);
    void bfs(std::size_t u);
    void printGraph();
    void transpose();
    std::size_t getCountNodesAtAGivenLevelBFS(std::size_t u, std::size_t level);
    std::size_t getCountNodesAtAGivenLevelDFS(std::size_t u, std::size_t level);
    std::size_t getCountOfAllPossiblePath(std::size_t u, std::size_t v);
    std::vector<std::size_t> topSort_Kahn();
    std::vector<std::size_t> topSort_dfs();
    bool is_cycled();
    void findSCC_Tarjan();
    void findSCC_Kosarajou();
    std::vector<int> findSSSP_dfs(std::size_t u);
    std::vector<std::size_t> findShortestPath_dfs(std::size_t u, std::size_t v);
    std::vector<int> findSSSP_dijkstra(std::size_t u);
private:
    std::size_t m_countOfvertices;
    std::vector<std::vector<int>> m_adjMatrix;
    bool m_is_directed;
    void dfs_helper(std::size_t u, std::vector<bool>& visited);
    void dfs_helper_countNodesAtAGivenLevel(std::size_t u, std::size_t level, std::size_t currentLevel, std::vector<bool>& visited, std::size_t& count);
    std::size_t dfs_helper_countOfAllPossiblePath(std::size_t u, std::size_t v, std::vector<bool>& visited);
    void dfs_helper_topSort(std::size_t u, std::vector<bool>& visited, std::vector<std::size_t>& topOrder);
    bool dfs_helper_isCycledUndirected(std::size_t u, std::vector<bool>& visited, std::size_t parent);
    bool dfs_helper_isCycledDirected(std::size_t u, std::vector<bool>& visited, std::vector<bool>& onStack);
    bool isCycledDirected();
    bool isCycledUndirected();
    void dfs_helper_Tarjan(std::size_t u, std::vector<int>& ids, std::vector<int>& lowLink, std::stack<std::size_t>& s, std::vector<bool>& onStack);
    void fillInOrder(std::size_t u, std::stack<std::size_t>& s, std::vector<bool>& visited);
    void dfs_helper_onTransposedGraph(std::size_t u, std::vector<bool>& visited);
    void dfs_helper_sssp(std::size_t u, std::vector<int>& dist);
    void dfs_helper_findShortestPath(std::size_t u, std::vector<int>& parent, std::vector<int>& dist);
};
#endif //GRAPH_h