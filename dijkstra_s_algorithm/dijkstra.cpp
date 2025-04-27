#include <iostream>
#include <vector>
#include <random>
#include <ctime>
#include <limits>

class Graph {
public:
    Graph(int size, double density, double distance_low, double distance_high)
        : m_size(size),
          m_density(density),
          m_distance_low(distance_low),
          m_distance_high(distance_high),
          m_graph(size, std::vector<double>(size, 0.0)),
          dist(size, std::numeric_limits<double>::max()),
          visited(size, false)
    {}

    void drawGraph();
    void printGraph() const;
    bool isConnected() const;
    void dijkstra() const;

private:
    int m_size;
    double m_density;
    double m_distance_low;
    double m_distance_high;
    std::vector<std::vector<double>> m_graph;
    mutable std::vector<double> dist;
    mutable std::vector<bool> visited;
};

void Graph::drawGraph() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(m_distance_low, m_distance_high);
    std::uniform_real_distribution<> prob(0.0, 1.0);

    for (int i = 0; i < m_size; ++i) {
        for (int j = i + 1; j < m_size; ++j) {
            if (prob(gen) < m_density) {
                double distance = dis(gen);
                m_graph[i][j] = distance;
                m_graph[j][i] = distance; // symmetric
            }
        }
    }
}

void Graph::printGraph() const {
    for (const auto& row : m_graph) {
        for (double val : row) {
            std::cout << val << "\t";
        }
        std::cout << "\n";
    }
}

bool Graph::isConnected() const {
    std::vector<bool> open(m_size, false);
    std::vector<bool> closed(m_size, false);
    int closed_size = 0;

    open[0] = true;

    while (closed_size < m_size) {
        int prev_closed_size = closed_size;
        for (int i = 0; i < m_size; ++i) {
            if (open[i] && !closed[i]) {
                closed[i] = true;
                ++closed_size;
                for (int j = 0; j < m_size; ++j) {
                    if (m_graph[i][j] > 0) {
                        open[j] = true;
                    }
                }
            }
        }
        if (closed_size == prev_closed_size) {
            return false; // no progress
        }
    }
    return true;
}

void Graph::dijkstra() const {
    dist.assign(m_size, std::numeric_limits<double>::max());
    visited.assign(m_size, false);

    dist[0] = 0.0;

    for (int count = 0; count < m_size - 1; ++count) {
        double min_dist = std::numeric_limits<double>::max();
        int min_index = -1;

        for (int i = 0; i < m_size; ++i) {
            if (!visited[i] && dist[i] < min_dist) {
                min_dist = dist[i];
                min_index = i;
            }
        }

        if (min_index == -1) {
            break; // Remaining nodes are unreachable
        }

        visited[min_index] = true;

        for (int i = 0; i < m_size; ++i) {
            if (!visited[i] && m_graph[min_index][i] > 0 &&
                dist[min_index] + m_graph[min_index][i] < dist[i]) {
                dist[i] = dist[min_index] + m_graph[min_index][i];
            }
        }
    }

    double sum = 0.0;
    std::cout << "Vertex\tDistance from Source\n";
    for (int i = 0; i < m_size; ++i) {
        std::cout << i << "\t" << dist[i] << "\n";
        sum += dist[i];
    }

    double average_dist = sum / m_size;
    std::cout << "Average distance: " << average_dist << "\n";
}

int main() {
    srand(static_cast<unsigned>(time(nullptr)));

    std::cout << "Draw a graph: 20% density, 50 nodes, distance range 1.0 to 10.0\n";
    Graph g1(50, 0.2, 1.0, 10.0);
    g1.drawGraph();
    g1.printGraph();

    if (g1.isConnected()) {
        std::cout << "The generated graph is connected.\n";
        g1.dijkstra();
    } else {
        std::cout << "The generated graph is not connected.\n";
    }

    std::cout << "\nDraw a graph: 40% density, 50 nodes, distance range 1.0 to 10.0\n";
    Graph g2(50, 0.4, 1.0, 10.0);
    g2.drawGraph();
    g2.printGraph();

    if (g2.isConnected()) {
        std::cout << "The generated graph is connected.\n";
        g2.dijkstra();
    } else {
        std::cout << "The generated graph is not connected.\n";
    }

    return 0;
}
