#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <iostream>
#include <cmath>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
 
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <omp.h>
#include <random>

#include "vector.h"
#include "scene.h"
#include "trianglemesh.h"


// The code for the boxmuller function was taken online
std::pair<double, double> boxMuller(double mean = 0.0, double stddev = 1.0) {
    static std::mt19937 rng(omp_get_thread_num());
    static std::uniform_real_distribution<double> uni01(0.0, 1.0);

    double u1 = uni01(rng);
    double u2 = uni01(rng);

    double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    double z1 = sqrt(-2.0 * log(u1)) * sin(2.0 * M_PI * u2);

    return {mean + z0 * stddev, mean + z1 * stddev};
}

int main() {
    int W = 1024;
    int H = 1024;
    double angle = 60 * M_PI / 180;

    int samplesPerPixel = 128;
    
    Vector Camera(0, 0, 110.0);
    Vector Light(-20.0, 40.0, 80.0);
    double LightIntensity = 6e9;

    Sphere* bottom = new Sphere(Vector(0, -2000.0, 0), 1980.0, Vector(0.3, 0.4, 0.7));
    Sphere* top = new Sphere(Vector(0, 2000.0, 0), 1880.0, Vector(0.2, 0.5, 0.9));
    Sphere* left = new Sphere(Vector(-2000.0, 0, 0), 1880.0, Vector(0.9, 0.2, 0.9));
    Sphere* right = new Sphere(Vector(2000.0, 0, 0), 1880.0, Vector(0.6, 0.5, 0.1));
    Sphere* back = new Sphere(Vector(0, 0, -2000.0), 1880.0, Vector(0.4, 0.8, 0.7));
    Sphere* front = new Sphere(Vector(0, 0, 2000.0), 1880.0, Vector(0.9, 0.4, 0.3));

    std::cerr << "Testing 1" << std::endl;

    TriangleMesh* cat = new TriangleMesh(Vector(1, 1, 1), NORMAL);
    cat->readOBJ("cadnav.com_model/Models_F0202A090/cat.obj");

    std::cerr << "Testing 2" << std::endl;

    for (auto& vertex : cat->vertices) {
        vertex = vertex * 1.2;     
        vertex = vertex + Vector(0, -20, 0); 
    } // I am not sure if this is a vary efficient way to scale and move the mash grid, but it works for now
    
    std::cerr << "Testing 2.1" << std::endl;

    cat->computeBoundingBox();

    std::cerr << "Testing 2.2" << std::endl;

    cat->buildBVH();

    std::cerr << "Testing 2.5" << std::endl;

    //Sphere hollowSphere1(Vector(20.0, 0, 0), 10.0, Vector(1.0, 1.0, 1.0), TRANSPARENT, 1.0, 1.5, false);
    //Sphere hollowSphere2(Vector(20.0, 0, 0), 9.9, Vector(1.0, 1.0, 1.0), TRANSPARENT, 1.0, 1.5, true);
    //Sphere sphere(Vector(20.0, 0, 0), 10.0, Vector(0.0, 1.0, 1.0), NORMAL);
    //Sphere mirrorsphere(Vector(-20.0, 0, 0), 10.0, Vector(1, 1, 1), MIRROR);
    //Sphere glassSphere(Vector(0, 0, 0), 10.0, Vector(1.0, 1.0, 1.0), TRANSPARENT, 1.0, 1.5);

    Scene scene(Light, LightIntensity);

    std::cerr << "Testing 3" << std::endl;

    scene.add(bottom);
    scene.add(top);
    scene.add(left);
    scene.add(right);
    scene.add(back);
    scene.add(front);

    std::cerr << "Testing 4" << std::endl;

    scene.add(cat);

    std::cerr << "Testing 5" << std::endl;

    //scene.add(sphere);
    //scene.add(mirrorsphere);
    //scene.add(glassSphere);

    std::vector<Vector> images2(W * H); 


    #pragma omp parallel for schedule(dynamic, 1)

    for (int i = 0; i < H; i++) {
        std::mt19937 rng(omp_get_thread_num());
        std::uniform_real_distribution<double> uni01(0.0, 1.0);

        for (int j = 0; j < W; j++) {

            Vector res(0, 0, 0);

            for (int s = 0; s < samplesPerPixel; ++s) {
                auto [x2, y2] = boxMuller(0.0, 0.5);
                double z = -(W / (2 * tan(angle / 2)));
                double x = (j - W / 2 + 0.5 + x2);
                double y = (H / 2 - i - 0.5 + y2);
            

                Vector dir(x, y, z);
                dir.normalize();
                Ray ray(Camera, dir);

                res = res + scene.getColor(ray, 5);
            }

            res = res / samplesPerPixel;

            images2[i * W + j] = res;
        }
    }

    // I saw this wasn't done the same in the source code, but it said to do the gamma correction at the very end so I am looping over again here.
    // although, it might not be very efficient

    std::vector<unsigned char> image(W * H * 3);

    for (int i = 0; i < W * H; i++) {
        for (int c = 0; c < 3; c++) {
            double base = std::max(0.0, images2[i][c]);
            double gamma_corrected = std::min(255.0, pow(base, 1.0 / 2.2));
            image[i * 3 + c] = static_cast<unsigned char>(gamma_corrected);
        }
    }

    stbi_write_png("image.png", W, H, 3, &image[0], 0);

    delete bottom;
    delete top;
    delete left;
    delete right;
    delete back;
    delete front;
    delete cat;

    return 0;
}