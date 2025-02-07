#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <cstdint>
#include <limits>
#include <sstream>
#include <zlib.h>  // Required for GZip compression
#include "dft_library.h"  // Placeholder for LLNL/INQ DFT package

#pragma pack(push, 1)
struct MRCHeader {
    int32_t nx, ny, nz;
    int32_t mode;
    int32_t nxstart, nystart, nzstart;
    int32_t mx, my, mz;
    float xlen, ylen, zlen;
    float alpha, beta, gamma;
    int32_t mapc, mapr, maps;
    float amin, amax, amean;
    int32_t ispg, nsymbt;
    char extra[100];
    float origin[3];
    char map[4];
    int32_t machst;
    float rms;
    int32_t nlabel;
    char labels[10][80];
};
#pragma pack(pop)

void normalize_density(std::vector<float>& density_data, float min_target = 0.0f, float max_target = 1.0f) {
    float min_val = *std::min_element(density_data.begin(), density_data.end());
    float max_val = *std::max_element(density_data.begin(), density_data.end());
    for (float& value : density_data) {
        value = min_target + (value - min_val) / (max_val - min_val) * (max_target - min_target);
    }
}

void extract_region(std::vector<float>& density_data, int nx, int ny, int nz,
                     int start_x, int end_x, int start_y, int end_y, int start_z, int end_z) {
    std::vector<float> extracted;
    for (int ix = start_x; ix < end_x; ix++) {
        for (int iy = start_y; iy < end_y; iy++) {
            for (int iz = start_z; iz < end_z; iz++) {
                extracted.push_back(density_data[ix + nx * (iy + ny * iz)]);
            }
        }
    }
    density_data = std::move(extracted);
}

void write_mrc(const std::string& filename, const DensityGrid& density, bool normalize = false) {
    auto basis = density.basis();
    MRCHeader header = {};
    header.nx = basis.sizes()[0];
    header.ny = basis.sizes()[1];
    header.nz = basis.sizes()[2];
    header.mode = 2;
    header.mx = header.nx;
    header.my = header.ny;
    header.mz = header.nz;
    header.xlen = basis.lattice_constants()[0];
    header.ylen = basis.lattice_constants()[1];
    header.zlen = basis.lattice_constants()[2];
    strncpy(header.map, "MAP ", 4);
    header.machst = 0x44;

    std::vector<float> density_data(header.nx * header.ny * header.nz);
    float min_val = std::numeric_limits<float>::max();
    float max_val = std::numeric_limits<float>::lowest();
    float sum = 0.0f;

    for (int ix = 0; ix < header.nx; ix++) {
        for (int iy = 0; iy < header.ny; iy++) {
            for (int iz = 0; iz < header.nz; iz++) {
                float value = density.cubic()[ix][iy][iz];
                density_data[ix + header.nx * (iy + header.ny * iz)] = value;
                min_val = std::min(min_val, value);
                max_val = std::max(max_val, value);
                sum += value;
            }
        }
    }

    if (normalize) normalize_density(density_data);
    
    header.amin = min_val;
    header.amax = max_val;
    header.amean = sum / density_data.size();

    std::ofstream outFile(filename, std::ios::binary);
    if (!outFile) {
        std::cerr << "Error: Cannot open file " << filename << " for writing.\n";
        return;
    }

    outFile.write(reinterpret_cast<const char*>(&header), sizeof(MRCHeader));
    outFile.write(reinterpret_cast<const char*>(density_data.data()), density_data.size() * sizeof(float));
    outFile.close();
    std::cout << "MRC file written: " << filename << std::endl;
}

void write_mrc_gzip(const std::string& filename, const DensityGrid& density, bool normalize = false) {
    gzFile file = gzopen(filename.c_str(), "wb");
    if (!file) {
        std::cerr << "Error: Cannot open file " << filename << " for writing.\n";
        return;
    }

    auto basis = density.basis();
    MRCHeader header = {};
    header.nx = basis.sizes()[0];
    header.ny = basis.sizes()[1];
    header.nz = basis.sizes()[2];
    header.mode = 2;
    header.mx = header.nx;
    header.my = header.ny;
    header.mz = header.nz;
    header.xlen = basis.lattice_constants()[0];
    header.ylen = basis.lattice_constants()[1];
    header.zlen = basis.lattice_constants()[2];
    strncpy(header.map, "MAP ", 4);
    header.machst = 0x44;

    std::vector<float> density_data(header.nx * header.ny * header.nz);
    for (int ix = 0; ix < header.nx; ix++) {
        for (int iy = 0; iy < header.ny; iy++) {
            for (int iz = 0; iz < header.nz; iz++) {
                density_data[ix + header.nx * (iy + header.ny * iz)] = density.cubic()[ix][iy][iz];
            }
        }
    }
    if (normalize) normalize_density(density_data);

    gzwrite(file, &header, sizeof(MRCHeader));
    gzwrite(file, density_data.data(), density_data.size() * sizeof(float));
    gzclose(file);
    std::cout << "GZ MRC file written: " << filename << std::endl;
}

int main() {
    auto density = electrons.density();
    write_mrc("electron_density.mrc", density, true);
    write_mrc_gzip("electron_density.mrc.gz", density, true);
    return 0;
}
