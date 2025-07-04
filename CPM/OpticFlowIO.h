#ifndef _OpticFlowIO_H
#define _OpticFlowIO_H

#include <cmath>
#include <cstring>
#include <cstdint>
#include <cstdio>
#include <algorithm>
#include <opencv2/opencv.hpp>

#define UNKNOWN_FLOW 1e10
#define UNKNOWN_FLOW_THRESH 1e9
#define TAG_FLOAT 202021.25f
#define TAG_STRING "PIEH"
#define NUM_BANDS 2
#define M_PI 3.14159265358979323846

typedef struct {
    double aee; // Average Endpoint Error
    double aae; // Average Angular Error
} FlowErr;

class OpticFlowIO {
public:
    template <class T>
    static bool unknown_flow(T u, T v) {
        return std::abs(u) > UNKNOWN_FLOW_THRESH || std::abs(v) > UNKNOWN_FLOW_THRESH || std::isnan(u) || std::isnan(v);
    }

    template <class T>
    static bool unknown_flow(T* f) {
        return unknown_flow(f[0], f[1]);
    }

    template <class T>
    static int ReadFlowFile(T* U, T* V, int* w, int* h, const char* filename) {
        if (!filename) return -1;

        FILE* stream = fopen(filename, "rb");
        if (!stream) return -1;

        float tag;
        int width, height;

        if (fread(&tag, sizeof(float), 1, stream) != 1 ||
            fread(&width, sizeof(int), 1, stream) != 1 ||
            fread(&height, sizeof(int), 1, stream) != 1 ||
            tag != TAG_FLOAT || width <= 0 || height <= 0) {
            fclose(stream);
            return -1;
        }

        *w = width;
        *h = height;

        for (int i = 0; i < width * height; ++i) {
            float tmp[2];
            if (fread(tmp, sizeof(float), 2, stream) != 2) {
                fclose(stream);
                return -1;
            }
            U[i] = tmp[0];
            V[i] = tmp[1];
        }

        fclose(stream);
        return 0;
    }

    template <class T>
    static int WriteFlowFile(T* U, T* V, int w, int h, const char* filename) {
        if (!filename) return -1;

        FILE* stream = fopen(filename, "wb");
        if (!stream) return -1;

        float tag = TAG_FLOAT;
        fwrite(&tag, sizeof(float), 1, stream);
        fwrite(&w, sizeof(int), 1, stream);
        fwrite(&h, sizeof(int), 1, stream);

        for (int i = 0; i < w * h; ++i) {
            float tmp[2] = { static_cast<float>(U[i]), static_cast<float>(V[i]) };
            fwrite(tmp, sizeof(float), 2, stream);
        }

        fclose(stream);
        return 0;
    }

    template <class T>
    static int ReadKittiFlowFile(T* U, T* V, int* w, int* h, const char* filename) {
        if (!filename || std::strstr(filename, ".png") == nullptr) return -1;

        cv::Mat img = cv::imread(filename, cv::IMREAD_UNCHANGED);
        if (img.empty() || img.type() != CV_16UC3) return -1;

        *w = img.cols;
        *h = img.rows;

        for (int y = 0; y < img.rows; ++y) {
            for (int x = 0; x < img.cols; ++x) {
                cv::Vec3w pix = img.at<cv::Vec3w>(y, x);
                if (pix[0] > 0) {
                    U[y * img.cols + x] = (pix[2] - 32768.0f) / 64.0f;
                    V[y * img.cols + x] = (pix[1] - 32768.0f) / 64.0f;
                } else {
                    U[y * img.cols + x] = UNKNOWN_FLOW;
                    V[y * img.cols + x] = UNKNOWN_FLOW;
                }
            }
        }

        return 0;
    }

    template <class T>
    static int WriteKittiFlowFile(T* U, T* V, int w, int h, const char* filename) {
        if (!filename || std::strstr(filename, ".png") == nullptr) return -1;

        cv::Mat img(h, w, CV_16UC3);

        for (int y = 0; y < h; ++y) {
            for (int x = 0; x < w; ++x) {
                int idx = y * w + x;
                if (!unknown_flow(U[idx], V[idx])) {
                    uint16_t u_enc = std::clamp(static_cast<int>(U[idx] * 64.0f + 32768.0f), 0, 65535);
                    uint16_t v_enc = std::clamp(static_cast<int>(V[idx] * 64.0f + 32768.0f), 0, 65535);
                    img.at<cv::Vec3w>(y, x) = cv::Vec3w(1, v_enc, u_enc);
                } else {
                    img.at<cv::Vec3w>(y, x) = cv::Vec3w(0, 0, 0);
                }
            }
        }

        std::vector<int> compression = { cv::IMWRITE_PNG_COMPRESSION, 1 };
        return cv::imwrite(filename, img, compression) ? 0 : -1;
    }

};

#endif // _OpticFlowIO_H
