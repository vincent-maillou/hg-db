#include "hypergraph_reorder/reorderer.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

using namespace hypergraph_reorder;

// Returns path without its extension, e.g. "/foo/bar.mtx" -> "/foo/bar"
static std::string stem(const std::string& path) {
    auto dot = path.find_last_of('.');
    return dot != std::string::npos ? path.substr(0, dot) : path;
}

static void write_permutation(const std::string& path, const std::vector<index_t>& perm) {
    std::ofstream f(path);
    if (!f) throw HypergraphReorderError("Cannot open file for writing: " + path);
    for (auto v : perm) f << v << "\n";
}

static void print_usage(const char* prog) {
    std::cout <<
        "Usage: " << prog << " <input_matrix> <k> [options]\n\n"
        "Arguments:\n"
        "  <input_matrix>   Input matrix file (.mtx, .graph, or binary CSR)\n"
        "  <k>              Number of diagonal blocks/parts\n\n"
        "Options:\n"
        "  -h, --help           Show this help message\n"
        "  --preset <name>      MT-KaHyPar preset: default, quality, deterministic, large_k\n"
        "  --threads <n>        Number of threads (default: auto)\n"
        "  --quiet              Suppress partitioner output\n"
        "  --output-perm        Write permutation to <input>.perm\n"
        "  --output-graph       Write METIS graph to <input>.graph\n";
}

static void print_statistics(const Statistics& s) {
    std::cout << "\n=== Statistics ===\n"
              << "Matrix:     " << s.n_rows << " x " << s.n_cols << ", " << s.nnz << " nnz\n"
              << "Graph:      " << s.n_vertices << " vertices, " << s.n_edges << " edges\n"
              << "Cliques:    " << s.n_cliques
              << " (max " << s.max_clique_size
              << ", avg " << std::fixed << std::setprecision(1) << s.avg_clique_size << ")\n"
              << "Hypergraph: " << s.n_hypernodes << " nodes, "
              << s.n_hyperedges << " nets, " << s.total_pins << " pins\n"
              << "Partition:  " << s.n_parts << " parts, separator "
              << s.separator_size << " (" << std::setprecision(2)
              << (s.separator_ratio * 100) << "%)\n"
              << "Ordering:   " << s.blocks_ordered << " blocks ("
              << s.blocks_failed << " failed)\n"
              << "Time:       " << s.time_total_ms << " ms\n";
}

static MtKahyparPreset parse_preset(const std::string& s) {
    if (s == "quality")       return MtKahyparPreset::QUALITY;
    if (s == "deterministic") return MtKahyparPreset::DETERMINISTIC;
    if (s == "large_k")       return MtKahyparPreset::LARGE_K;
    if (s != "default")
        std::cerr << "Warning: unknown preset '" << s << "', using default\n";
    return MtKahyparPreset::DEFAULT;
}

int main(int argc, char** argv) {
    if (argc < 2 || std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help") {
        print_usage(argv[0]);
        return argc < 2 ? 1 : 0;
    }
    if (argc < 3) {
        std::cerr << "Error: expected <input_matrix> and <k>\n\n";
        print_usage(argv[0]);
        return 1;
    }

    const std::string input_path = argv[1];
    const int n_parts = std::atoi(argv[2]);
    if (n_parts < 2) { std::cerr << "Error: k must be >= 2\n"; return 1; }

    const std::string base = stem(input_path);
    const std::string output_path = base + "_reordered.mtx";

    MtKahyparPreset preset = MtKahyparPreset::DEFAULT;
    int  num_threads  = 0;
    bool quiet        = false;
    bool output_perm  = false;
    bool output_graph = false;

    for (int i = 3; i < argc; ++i) {
        std::string arg = argv[i];
        if      (arg == "--preset"  && i+1 < argc) preset      = parse_preset(argv[++i]);
        else if (arg == "--threads" && i+1 < argc) num_threads = std::atoi(argv[++i]);
        else if (arg == "--quiet")                 quiet        = true;
        else if (arg == "--output-perm")           output_perm  = true;
        else if (arg == "--output-graph")          output_graph = true;
        else if (arg == "-h" || arg == "--help") { print_usage(argv[0]); return 0; }
        else { std::cerr << "Error: unknown option '" << arg << "'\n\n"; print_usage(argv[0]); return 1; }
    }

    SymmetricDBReorderer::Options opts;
    opts.n_parts                    = n_parts;
    opts.imbalance                  = 0.03;
    opts.preset                     = preset;
    opts.seed                       = -1;
    opts.ordering_method            = OrderingMethod::AMD;
    opts.use_openmp                 = true;
    opts.num_threads                = num_threads;
    opts.suppress_partitioner_output = quiet;

    try {
        CSRMatrix original_matrix;
        if (output_graph)
            original_matrix = read_matrix(input_path, MatrixFormat::AUTO);

        SymmetricDBReorderer reorderer(opts);
        auto result = reorderer.reorder_from_file(input_path, MatrixFormat::AUTO);

        std::cout << "Saving reordered matrix to " << output_path << "\n";
        write_matrix(output_path, result.reordered_matrix);

        if (output_perm || output_graph) {
            std::string perm_path = base + ".perm";
            std::cout << "Saving permutation to " << perm_path << "\n";
            write_permutation(perm_path, result.permutation);
        }

        if (output_graph) {
            std::string graph_path = base + ".graph";
            std::cout << "Saving METIS graph to " << graph_path << "\n";
            metis::write(graph_path, original_matrix);
        }

        print_statistics(result.stats);
        std::cout << "\nDone.\n";
        return 0;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
