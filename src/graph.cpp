#include "../include/graph.h"

Graph::Graph(std::size_t countOfVertices, bool is_directed) : m_countOfvertices(countOfVertices), m_is_directed(is_directed)
{
    this->m_adjMatrix.resize(this->m_countOfvertices, std::vector<int>(this->m_countOfvertices, std::numeric_limits<int>::min()));
}

void Graph::add_edge(std::size_t u, std::size_t v, int weight)
{
    if(u < this->m_countOfvertices && v < this->m_countOfvertices)
    {
        this->m_adjMatrix[u][v] = weight;
        if(!this->m_is_directed)
        {
            this->m_adjMatrix[v][u] = weight;
        }
    }
}

void Graph::remove_edge(std::size_t u, std::size_t v)
{
    if(u < this->m_countOfvertices && v < this->m_countOfvertices)
    {
        this->m_adjMatrix[u][v] = std::numeric_limits<int>::min();
        if(!this->m_is_directed)
        {
            this->m_adjMatrix[v][u] = std::numeric_limits<int>::min();
        }
    }
}

void Graph::dfs(std::size_t u)
{
    std::vector<bool> visited(this->m_countOfvertices, false);
    dfs_helper(u, visited);
    std::cout << std::endl;
}

void Graph::dfs_helper(std::size_t u, std::vector<bool>& visited)
{
    visited[u] = true;
    std::cout << u << " ";
    for(int v = 0; v < this->m_countOfvertices; ++v)
    {
        if(visited[v] && this->m_adjMatrix[u][v] != std::numeric_limits<int>::min())
        {
            dfs_helper(v, visited);
        }
    }
}

void Graph::bfs(std::size_t u)
{
    std::vector<bool> visited(this->m_countOfvertices, false);
    std::queue<std::size_t> q;
    q.push(u);
    visited[u] = true;
    while(!q.empty())
    {
        int current = q.front();
        std::cout << current << " ";
        q.pop();

        for(int v = 0; v < this->m_countOfvertices; ++v)
        {
            if(!visited[v] && this->m_adjMatrix[u][v] != std::numeric_limits<int>::min())
            {
                q.push(v);
                visited[v] = true;
            }
        }

    }
    std::cout << std::endl;
}

void Graph::printGraph()
{
    for(int u = 0; u < this->m_countOfvertices; ++u)
    {
        for(int v = 0; v < this->m_countOfvertices; ++v)
        {
            if(this->m_adjMatrix[u][v] == std::numeric_limits<int>::min())
            {
                std::cout << "N" << " ";
            } else 
            {
                std::cout << this->m_adjMatrix[u][v] << " ";
            }
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void Graph::transpose()
{
    std::vector<std::vector<int>> transposedGraph(this->m_countOfvertices, std::vector<int>(this->m_countOfvertices, std::numeric_limits<int>::min()));
    for(int u = 0; u < this->m_countOfvertices; ++u)
    {
        for(int v = 0; v < this->m_countOfvertices; ++v)
        {
            transposedGraph[v][u] = this->m_adjMatrix[u][v];
        }
    }
    this->m_adjMatrix = transposedGraph;
}

std::size_t Graph::getCountNodesAtAGivenLevelBFS(std::size_t u, std::size_t level)
{
    std::vector<bool> visited(this->m_countOfvertices, false);
    std::size_t count = 0;
    std::vector<int> levels(this->m_countOfvertices, 0);
    std::queue<std::size_t> q;
    q.push(u);
    visited[u] = true;
    while(!q.empty())
    {
        std::size_t current = q.front();
        q.pop();
        for(int v = 0; v < this->m_countOfvertices; ++v)
        {
            if(!visited[current] && this->m_adjMatrix[current][v] != std::numeric_limits<int>::min())
            {
                q.push(v);
                visited[v] = true;
                levels[v] = levels[current] + 1;
            }
        }
    }

    std::size_t levelsSize = levels.size();
    for(int i = 0; i < levelsSize; ++i)
    {
        if(i == level) ++count;
    }
    return count;
}

std::size_t Graph::getCountNodesAtAGivenLevelDFS(std::size_t u, std::size_t level)
{
    std::vector<bool> visited(this->m_countOfvertices, false);
    std::size_t count = 0;
    dfs_helper_countNodesAtAGivenLevel(u, level, 0, visited, count);
    return count;
}

void Graph::dfs_helper_countNodesAtAGivenLevel(std::size_t u, std::size_t level, std::size_t currentLevel, std::vector<bool>& visited, std::size_t& count)
{
    if(currentLevel == level) ++count;

    for(int v = 0; v < this->m_countOfvertices; ++v)
    {
        if(!visited[v] && this->m_adjMatrix[u][v] != std::numeric_limits<int>::min())
        {
            dfs_helper_countNodesAtAGivenLevel(v, level, currentLevel + 1, visited, count);
        }
    }
}

std::size_t Graph::dfs_helper_countOfAllPossiblePath(std::size_t u, std::size_t v, std::vector<bool>& visited)
{
    if(u == v) return 1;

    visited[u] = true;
    std::size_t count = 0;
    for(std::size_t i = 0; i < this->m_countOfvertices; ++i)
    {
        if(!visited[i] && this->m_adjMatrix[u][i] != std::numeric_limits<int>::min())
        {
            count += dfs_helper_countOfAllPossiblePath(i, v, visited);
        }
    }
    visited[u] = false;
    return count;
}

std::size_t Graph::getCountOfAllPossiblePath(std::size_t u, std::size_t v)
{
    std::vector<bool> visited(this->m_countOfvertices, false);
    return dfs_helper_countOfAllPossiblePath(u, v, visited);
}

std::vector<std::size_t> Graph::topSort_Kahn()
{
    std::vector<std::size_t> inDegree(this->m_countOfvertices, 0);
    for(int u = 0; u < this->m_countOfvertices; ++u)
    {
        for(int v = 0; v < this->m_countOfvertices; ++v)
        {
            if(this->m_adjMatrix[u][v] != std::numeric_limits<int>::min())
            {
                ++inDegree[v];
            }
        }
    }

    std::queue<std::size_t> q;
    std::vector<std::size_t> topOrder;
    for(int i = 0; i < this->m_countOfvertices; ++i)
    {
        if(inDegree[i] == 0)
        {
            q.push(i);
        }
    }

    while (!q.empty())
    {
        std::size_t current = q.front();
        topOrder.push_back(current);
        q.pop();
        for(int v = 0; v < this->m_countOfvertices; ++v)
        {
            if(this->m_adjMatrix[current][v] != std::numeric_limits<int>::min())
            {
                --inDegree[v];
                if(inDegree[v] == 0)
                {
                    q.push(v);
                }
            }
        }
    }
    if(this->m_countOfvertices != topOrder.size()) throw std::logic_error("Graph is Cycled. TopSort is impossible");
    return topOrder;
    
}

void Graph::dfs_helper_topSort(std::size_t u, std::vector<bool>& visited, std::vector<std::size_t>& topOrder)
{
    visited[u] = true;
    for(int v = 0; v < this->m_countOfvertices; ++v)
    {
        if(!visited[v] && this->m_adjMatrix[u][v] != std::numeric_limits<int>::min())
        {
            dfs_helper_topSort(v, visited, topOrder);
        }
    }
    topOrder.push_back(u);
}

std::vector<std::size_t> Graph::topSort_dfs()
{
    std::vector<bool> visited(this->m_countOfvertices, false);
    std::vector<std::size_t> topOrder;
    for(int u = 0; u < this->m_countOfvertices; ++u)
    {
        if(!visited[u]) 
        {
            dfs_helper_topSort(u, visited, topOrder);
        }
    }
    if(topOrder.size() != this->m_countOfvertices) throw std::logic_error("Graph is Cycled. TopSort is impossible");
    std::reverse(topOrder.begin(), topOrder.end());
    return topOrder;
}

bool Graph::dfs_helper_isCycledUndirected(std::size_t u, std::vector<bool>& visited, std::size_t parent)
{
    visited[u] = true;
    for(int v = 0; v < this->m_countOfvertices; ++v)
    {
        if(v != parent && this->m_adjMatrix[u][v] != std::numeric_limits<int>::min())
        {
            if(visited[v]) return true;
            if(dfs_helper_isCycledUndirected(v, visited, parent)) return true;
        }
    }
    return false;
}

bool Graph::dfs_helper_isCycledDirected(std::size_t u, std::vector<bool>& visited, std::vector<bool>& onStack)
{
    visited[u] = true;
    onStack[u] = true;

    for(int v = 0; v < this->m_countOfvertices; ++v)
    {
        if(!visited[v] && this->m_adjMatrix[u][v] != std::numeric_limits<int>::min())
        {
            if(dfs_helper_isCycledDirected(v, visited, onStack)) return true;
        } else if(onStack[v])
        {
            return true;
        }
    }
    onStack[u] = false;
    return false;
}

bool Graph::isCycledDirected()
{
    std::vector<bool> visited(this->m_countOfvertices, false);
    std::vector<bool> onStack(this->m_countOfvertices, false);
    for(int u = 0; u < this->m_countOfvertices; ++u)
    {
        if(!visited[u])
        {
            if(dfs_helper_isCycledDirected(u, visited, onStack)) return true;
        }
    }
    return false;
}

bool Graph::isCycledUndirected()
{
    std::vector<bool> visited(this->m_countOfvertices, false);
    for(int u = 0; u < this->m_countOfvertices; ++u)
    {
        if(!visited[u])
        {
            if(dfs_helper_isCycledUndirected(u, visited, -1)) return true;
        }
    }
    return false;
}

bool Graph::is_cycled()
{
    return this->m_is_directed ? isCycledDirected() : isCycledUndirected();
}

void Graph::dfs_helper_Tarjan(std::size_t u, std::vector<int>& ids, std::vector<int>& lowLink, std::stack<std::size_t>& s, std::vector<bool>& onStack)
{
    static int id = 0;
    ids[u] = lowLink[u] = id++;
    onStack[u] = true;
    s.push(u);
    for(int v = 0; v < this->m_countOfvertices; ++v)
    {
        if(this->m_adjMatrix[u][v] != std::numeric_limits<int>::min())
        {
            if(ids[v] == -1)
            {
                dfs_helper_Tarjan(v, ids, lowLink, s, onStack);
            }
            if(onStack[v])
            {
                lowLink[u] = std::min(lowLink[u], lowLink[v]);
            }
        }
    }
    if(lowLink[u] == ids[u])
    {
        std::cout << "SCC: ";
        while(true)
        {
            std::size_t current = s.top();
            s.pop();
            std::cout << current << " ";
            onStack[current] = false;
            if(current == u) break;
        }
        std::cout << std::endl;
    }
}

void Graph::findSCC_Tarjan()
{
    std::vector<int> ids(this->m_countOfvertices, -1);
    std::vector<int> lowLink(this->m_countOfvertices, -1);
    std::stack<std::size_t> s;
    std::vector<bool> onStack(this->m_countOfvertices, false);
    for(int i = 0; i < this->m_countOfvertices; ++i)
    {
        if(ids[i] == -1)
        {
            dfs_helper_Tarjan(i, ids, lowLink, s, onStack);
        }
    }
}

void Graph::fillInOrder(std::size_t u, std::stack<std::size_t>& s, std::vector<bool>& visited)
{
    visited[u] = true;
    for(int v = 0; v < this->m_countOfvertices; ++v)
    {
        if(visited[v] && this->m_adjMatrix[u][v] != std::numeric_limits<int>::min())
        {
            fillInOrder(v, s, visited);
        }
    }
    s.push(u);
}

void Graph::dfs_helper_onTransposedGraph(std::size_t u, std::vector<bool>& visited)
{
    visited[u] = true;
    std::cout << u << " ";
    for(int v = 0; v < this->m_countOfvertices; ++v)
    {
        if(!visited[v] && this->m_adjMatrix[u][v] != std::numeric_limits<int>::min())
        {
            dfs_helper_onTransposedGraph(v, visited);
        }
    }
}

void Graph::findSCC_Kosarajou()
{
    std::vector<bool> visited(this->m_countOfvertices, false);
    std::stack<std::size_t> s;
    for(int u = 0; u < this->m_countOfvertices; ++u)
    {
        if(!visited[u])
        {
            fillInOrder(u, s, visited);
        }
    }

    transpose();
    std::fill(visited.begin(), visited.end(), false);
    while(!s.empty())
    {

        std::size_t current = s.top();
        s.pop();
        if(!visited[current])
        {
            std::cout << "SCC: ";
            dfs_helper_onTransposedGraph(current, visited);
        }
    }
    std::cout << std::endl;
}

void Graph::dfs_helper_sssp(std::size_t u, std::vector<int>& dist)
{
    for(int v = 0; v < this->m_countOfvertices; ++v)
    {
        if(this->m_adjMatrix[u][v] != std::numeric_limits<int>::min())
        {
            int newDist = dist[u] + this->m_adjMatrix[u][v];
            if(newDist < dist[v])
            {
                dist[v] = newDist; 
                dfs_helper_sssp(v, dist);
            }
        }
    }
}

std::vector<int> Graph::findSSSP_dfs(std::size_t u)
{
    std::vector<bool> visited(this->m_countOfvertices, false);
    std::vector<int> dist(this->m_countOfvertices, std::numeric_limits<int>::max());
    dist[u] = 0;
    dfs_helper_sssp(u, dist);

    return dist;
}


void Graph::dfs_helper_findShortestPath(std::size_t u, std::vector<int>& parent, std::vector<int>& dist)
{
    for(int v = 0; v < this->m_countOfvertices; ++v)
    {
        if(this->m_adjMatrix[u][v] != std::numeric_limits<int>::max())
        {
            int newDist = dist[u] + this->m_adjMatrix[u][v];
            if(newDist < dist[v])
            {
                parent[u] = v;
                dist[v] = newDist;
                dfs_helper_findShortestPath(u, parent, dist);
            }
        }
    }
}

std::vector<std::size_t> Graph::findShortestPath_dfs(std::size_t u, std::size_t v)
{
    std::vector<int> dist(this->m_countOfvertices, std::numeric_limits<int>::max());
    std::vector<int> parent(this->m_countOfvertices, -1);
    std::vector<std::size_t> path;
    dist[u] = 0;
    dfs_helper_findShortestPath(u, parent, dist);
    if(dist[v] == std::numeric_limits<int>::max())
    {
        return {};
    }

    for(auto at = v; at != -1; at = parent[at])
    {
        path.push_back(at);
    }

    std::reverse(path.begin(), path.end());
    return path; 
}

std::vector<int> Graph::findSSSP_dijkstra(std::size_t u)
{
    std::vector<int> dist(this->m_countOfvertices, std::numeric_limits<int>::max());
    dist[u] = 0;
    std::priority_queue<std::pair<int, std::size_t>, std::vector<std::pair<int, std::size_t>>, std::greater<std::pair<int, std::size_t>>> pq;
    pq.push({0, u});
    while(!pq.empty())
    {
        auto topElement = pq.top();
        std::size_t current = topElement.second;
        int currentDist = topElement.first;
        pq.pop();
        if(currentDist > dist[current]) continue;

        for(std::size_t v = 0; v < this->m_countOfvertices; ++v)
        {
            if(this->m_adjMatrix[current][v] != std::numeric_limits<int>::min())
            {
                int newDist = currentDist + this->m_adjMatrix[current][v];
                if(newDist < dist[v])
                {
                    dist[v] = newDist;
                    pq.push({newDist, v});
                }
            }
        }
    }
    return dist;
}