#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <random>
#include <limits>
#include <cmath>

// Функция для чтения данных из CSV файла
std::vector<std::vector<double>> readDataFromCSV(const std::string& filename) {
    std::vector<std::vector<double>> data;
    std::ifstream file(filename);
    std::string line;
    while (std::getline(file, line)) {
        std::stringstream linestream(line);
        std::string cell;
        std::vector<double> parsedRow;
        while (std::getline(linestream, cell, ',')) {
            parsedRow.push_back(std::stod(cell));
        }
        data.push_back(parsedRow);
    }
    return data;
}

// Функция для записи кластеров в CSV файл
void writeClustersToCSV(const std::string& filename, const std::vector<std::vector<double>>& clusters) {
    std::ofstream file(filename);
    for (const auto& cluster : clusters) {
        for (size_t i = 0; i < cluster.size(); ++i) {
            file << cluster[i];
            if (i < cluster.size() - 1) file << ",";
        }
        file << "\n";
    }
}

// Функция для вычисления Евклидова расстояния между двумя точками
double euclideanDistance(const std::vector<double>& a, const std::vector<double>& b) {
    double sum = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        sum += (a[i] - b[i]) * (a[i] - b[i]);
    }
    return std::sqrt(sum);
}

std::vector<std::vector<double>> kMeansClustering(const std::vector<std::vector<double>>& data, int k);

void interpretClusters(const std::vector<std::vector<double>>& clusters);

int main() {
    std::string inputFilename = "data.csv";
    std::string outputFilename = "clusters.csv";
    int k = 0;

    std::cout << "Введите количество кластеров k: ";
    std::cin >> k;

    try {
        std::vector<std::vector<double>> data = readDataFromCSV(inputFilename);

        std::vector<std::vector<double>> clusters = kMeansClustering(data, k);
        interpretClusters(clusters);

        writeClustersToCSV(outputFilename, clusters);

        std::cout << "Кластеризация завершена. Результаты записаны в файл " << outputFilename << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Произошла ошибка: " << e.what() << std::endl;
    }

    return 0;
}

void kMeansClustering(std::vector<std::vector<double> >& data, int k, std::vector<int>& assignments) {
    size_t dimensions = data[0].size();
    std::vector<std::vector<double> > centroids(k, std::vector<double>(dimensions, 0));
    std::vector<std::vector<double> > oldCentroids(k, std::vector<double>(dimensions, 0));
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, data.size() - 1);
    for (int i = 0; i < k; ++i) {
        centroids[i] = data[dis(gen)];
    }
    bool centroidsChanged = true;
    while (centroidsChanged) {
        for (size_t i = 0; i < data.size(); ++i) {
            double minDistance = std::numeric_limits<double>::max();
            int clusterIndex = 0;
            for (int j = 0; j < k; ++j) {
                double distance = euclideanDistance(data[i], centroids[j]);
                if (distance < minDistance) {
                    minDistance = distance;
                    clusterIndex = j;
                }
            }
            assignments[i] = clusterIndex;
        }
        oldCentroids = centroids;
        std::vector<int> counts(k, 0);
        centroids.assign(k, std::vector<double>(dimensions, 0));
        for (size_t i = 0; i < data.size(); ++i) {
            int clusterIndex = assignments[i];
            std::transform(centroids[clusterIndex].begin(), centroids[clusterIndex].end(), data[i].begin(), centroids[clusterIndex].begin(), std::plus<double>());
            counts[clusterIndex]++;
        }
        for (int i = 0; i < k; ++i) {
            if (counts[i] != 0) {
                std::transform(centroids[i].begin(), centroids[i].end(), centroids[i].begin(),
            }
        }
        centroidsChanged = false;
        for (int i = 0; i < k; ++i) {
            if (euclideanDistance(centroids[i], oldCentroids[i]) > 0.001) {
                centroidsChanged = true;
                break;
            }
        }
    }
}

void interpretClusters(const std::vector<std::vector<double>>& clusters) {
    for (size_t i = 0; i < clusters.size(); ++i) {
        std::cout << "Кластер " << i + 1 << ": \n";
        std::vector<double> centroid(clusters[i].size(), 0);
        for (const auto& point : clusters) {
            for (size_t j = 0; j < point.size(); ++j) {
                centroid[j] += point[j];
            }
        }
        for (double& val : centroid) {
            val /= clusters[i].size();
        }

        std::cout << "Центроид: ";
        for (const auto& val : centroid) {
            std::cout << val << " ";
        }
        std::cout << "\n";

        double maxDistance = 0.0;
        for (const auto& point : clusters) {
            double distance = euclideanDistance(point, centroid);
            maxDistance = std::max(maxDistance, distance);
        }

        std::cout << "Максимальное расстояние до центроида: " << maxDistance << "\n";
        std::cout << "Количество точек в кластере: " << clusters[i].size() << "\n\n";
    }
}