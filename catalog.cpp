#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include<algorithm>

int decode_value(int base, const std::string& value) {
    return std::stoi(value, nullptr, base);  
}


std::vector<double> gaussian_elimination(std::vector<std::vector<double>> A, std::vector<double> b) {
    int n = A.size();
   
    for (int i = 0; i < n; ++i) {
     
        for (int j = i + 1; j < n; ++j) {
            if (fabs(A[j][i]) > fabs(A[i][i])) {
                std::swap(A[i], A[j]);
                std::swap(b[i], b[j]);
            }
        }

        for (int j = i + 1; j < n; ++j) {
            double factor = A[j][i] / A[i][i];
            for (int k = i; k < n; ++k) {
                A[j][k] -= factor * A[i][k];
            }
            b[j] -= factor * b[i];
        }
    }

    
    std::vector<double> x(n);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = b[i];
        for (int j = i + 1; j < n; ++j) {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }

    return x;
}


double find_constant_term(const std::map<std::string, std::map<std::string, std::string>>& json_input) {
    
    int n = std::stoi(json_input.at("keys").at("n"));
    int k = std::stoi(json_input.at("keys").at("k"));

    std::vector<std::pair<int, int>> points; 
    for (const auto& entry : json_input) {
        if (entry.first != "keys") {  
            int x = std::stoi(entry.first);  
            int base = std::stoi(entry.second.at("base"));
            std::string y_encoded = entry.second.at("value");
            int y = decode_value(base, y_encoded);
            points.push_back({x, y});
        }
    }

    std::sort(points.begin(), points.end());

    std::vector<std::vector<double>> A(k, std::vector<double>(k));  
    std::vector<double> b(k);  
    for (int i = 0; i < k; ++i) {
        int x = points[i].first;
        int y = points[i].second;
        A[i][0] = x * x;  
        A[i][1] = x;      
        A[i][2] = 1;      
        b[i] = y;         
    }

    std::vector<double> solution = gaussian_elimination(A, b);

    return solution[2];  
}

int main() {

    std::map<std::string, std::map<std::string, std::string>> json_input = {
        {"keys", {{"n", "4"}, {"k", "3"}}},
        {"1", {{"base", "10"}, {"value", "4"}}},
        {"2", {{"base", "2"}, {"value", "111"}}},
        {"3", {{"base", "10"}, {"value", "12"}}},
        {"6", {{"base", "4"}, {"value", "213"}}}
    };

    double constant_term = find_constant_term(json_input);
    std::cout << "The constant term (c) is: " << constant_term << std::endl;

    return 0;
}
