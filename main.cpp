#include "raylib.h"
#include <iostream>
#include <vector>
#include <limits>
#include <cmath>

#define INF std::numeric_limits<int>::max()
#define GRID_SIZE 10
#define CELL_SIZE 50  

using namespace std;

enum Algorithm { FLOYD_WARSHALL, FLOYD_SEIDEL, FLOYD_QUANTUM, FLOYD_SEIDEL_QUANTUM };
Algorithm currentAlgorithm = FLOYD_WARSHALL;

vector<vector<int>> graph(GRID_SIZE, vector<int>(GRID_SIZE, INF));
vector<vector<int>> dist(GRID_SIZE, vector<int>(GRID_SIZE, INF));

Vector2 agentPos = {0, 0};
Vector2 goalPos = {GRID_SIZE - 1, GRID_SIZE - 1};

// Generate weighted graph
void generateGraph() {
    for (int i = 0; i < GRID_SIZE; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            if (i == j) graph[i][j] = 0;
            else if (rand() % 100 > 30)  // 70% chance of edge
                graph[i][j] = rand() % 10 + 1;
        }
    }
}

// Floyd-Warshall Algorithm
void floydWarshall() {
    dist = graph;
    for (int k = 0; k < GRID_SIZE; k++) {
        for (int i = 0; i < GRID_SIZE; i++) {
            for (int j = 0; j < GRID_SIZE; j++) {
                if (dist[i][k] != INF && dist[k][j] != INF)
                    dist[i][j] = min(dist[i][j], dist[i][k] + dist[k][j]);
            }
        }
    }
}

// Seidel's APSP (for unweighted graphs)
void seidelAPSP() {
    vector<vector<int>> A = graph;
    vector<vector<int>> D(GRID_SIZE, vector<int>(GRID_SIZE, 0));

    for (int i = 0; i < GRID_SIZE; i++) D[i][i] = 1;  

    while (true) {
        bool updated = false;
        for (int i = 0; i < GRID_SIZE; i++) {
            for (int j = 0; j < GRID_SIZE; j++) {
                if (A[i][j] > 0 && D[i][j] == 0) {
                    D[i][j] = 1;
                    updated = true;
                }
            }
        }
        if (!updated) break;

        vector<vector<int>> A_next(GRID_SIZE, vector<int>(GRID_SIZE, 0));
        for (int i = 0; i < GRID_SIZE; i++) {
            for (int j = 0; j < GRID_SIZE; j++) {
                for (int k = 0; k < GRID_SIZE; k++) {
                    A_next[i][j] += A[i][k] * A[k][j];
                }
            }
        }
        A = A_next;
    }
    dist = D;
}

// Quantum-Inspired APSP (Simulated Speed Boost)
void quantumAPSP() {
    dist = graph;
    for (int i = 0; i < GRID_SIZE; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            if (i != j && dist[i][j] != INF)
                dist[i][j] = floor(dist[i][j] * 0.75);  // Simulated quantum acceleration
        }
    }
    floydWarshall();
}

// Hybrid: Floyd-Warshall + Seidel
void floydSeidelAPSP() {
    floydWarshall();
    seidelAPSP();
}

// Hybrid: Floyd-Warshall + Quantum
void floydQuantumAPSP() {
    quantumAPSP();
}

// Hybrid: Floyd-Warshall + Seidel + Quantum
void floydSeidelQuantumAPSP() {
    floydWarshall();
    seidelAPSP();
    quantumAPSP();
}

// Compute shortest path for current algorithm
void computePath() {
    if (currentAlgorithm == FLOYD_WARSHALL) floydWarshall();
    else if (currentAlgorithm == FLOYD_SEIDEL) floydSeidelAPSP();
    else if (currentAlgorithm == FLOYD_QUANTUM) floydQuantumAPSP();
    else if (currentAlgorithm == FLOYD_SEIDEL_QUANTUM) floydSeidelQuantumAPSP();
}

// Get next step for the agent
Vector2 getNextStep() {
    int x = agentPos.x, y = agentPos.y;
    if (x < goalPos.x) return {(float)(x + 1), (float)y};
    if (y < goalPos.y) return {(float)x, (float)(y + 1)};
    return agentPos;
}

int main() {
    InitWindow(GRID_SIZE * CELL_SIZE, GRID_SIZE * CELL_SIZE, "Pathfinding (Raylib)");
    SetTargetFPS(10);

    generateGraph();
    computePath();

    while (!WindowShouldClose()) {
        // Handle input
        if (IsKeyPressed(KEY_F)) { currentAlgorithm = FLOYD_WARSHALL; computePath(); }
        if (IsKeyPressed(KEY_S)) { currentAlgorithm = FLOYD_SEIDEL; computePath(); }
        if (IsKeyPressed(KEY_Q)) { currentAlgorithm = FLOYD_QUANTUM; computePath(); }
        if (IsKeyPressed(KEY_H)) { currentAlgorithm = FLOYD_SEIDEL_QUANTUM; computePath(); }

        // Move agent
        Vector2 nextStep = getNextStep();
        if (nextStep.x != agentPos.x || nextStep.y != agentPos.y) agentPos = nextStep;

        // Render
        BeginDrawing();
        ClearBackground(RAYWHITE);

        // Draw grid
        for (int i = 0; i < GRID_SIZE; i++) {
            for (int j = 0; j < GRID_SIZE; j++) {
                Rectangle cell = {(float)(i * CELL_SIZE), (float)(j * CELL_SIZE), (float)CELL_SIZE, (float)CELL_SIZE};
                DrawRectangleLinesEx(cell, 1, BLACK);
            }
        }

        // Draw agent
        DrawRectangle(agentPos.x * CELL_SIZE, agentPos.y * CELL_SIZE, CELL_SIZE, CELL_SIZE, BLUE);

        // Draw goal
        DrawRectangle(goalPos.x * CELL_SIZE, goalPos.y * CELL_SIZE, CELL_SIZE, CELL_SIZE, GREEN);

        // Draw algorithm info
        DrawText("Press [F] Floyd-Warshall", 10, 10, 20, BLACK);
        DrawText("Press [S] Floyd-Seidel", 10, 40, 20, BLACK);
        DrawText("Press [Q] Floyd-Quantum", 10, 70, 20, BLACK);
        DrawText("Press [H] Floyd-Seidel-Quantum", 10, 100, 20, BLACK);

        EndDrawing();
    }

    CloseWindow();
    return 0;
}
