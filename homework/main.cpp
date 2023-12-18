#include <QString>
#include <QtGui>
#include <QApplication>
#include <QImage>
#include <QLabel>
#include <QTextStream>

#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>

const double INTENSITY_MIN = 0.0;
const double INTENSITY_MAX = 255.0;
const int KERNEL_SIZE = 9;

std::vector<std::vector<double>> generate_kernel(int, double);
QImage apply_filter(const QImage&, const std::vector<std::vector<double>>&);

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "no file name provided" << std::endl;
        std::cerr << "Usage: <program> <path/to/grayscale_image>" << std::endl;
        return 1;
    }

    QString file_name(QString::fromLocal8Bit(argv[1]));
    QApplication app(argc, argv);
    QImage image;
    QLabel label;

    if (!image.load(file_name)) {
        QTextStream(stderr) << "Failed to load image: " << file_name << Qt::endl;
        return 1;
    }

    if (image.format() != QImage::Format_Grayscale8) {
        QTextStream(stderr) << "Image \"" << file_name << "\" not in Grayscale format" << Qt::endl;
        return 1;
    }

    //coefficient of how smooth the Gaussian curve to be
    // determines how blurry (farther from 0) or sharp (closer to 0) the image becomes
    double sigma = 3.0;
    auto kernel = generate_kernel(KERNEL_SIZE, sigma);

    QTextStream(stdout) << "Kernel: " << Qt::endl;
    for (const auto& elem : kernel) { 
        for (const auto& inner_elem : elem) {
            QTextStream(stdout) << inner_elem << "\t"; 
        }
        QTextStream(stdout) << "" << Qt::endl; 
    }

    QImage out_image = apply_filter(image, kernel);

    label.setPixmap(QPixmap::fromImage(out_image));
    label.show();

    return app.exec();
}

std::vector<std::vector<double>> generate_kernel(int size, double sigma) {
    std::vector<std::vector<double>> kernel(size, std::vector<double>(size));
    double sum = 0.0;
    int half_size = size / 2;

    for (int row_idx = -half_size; row_idx <= half_size; ++row_idx) {
        for (int col_idx = -half_size; col_idx <= half_size; ++col_idx) {
            kernel[row_idx + half_size][col_idx + half_size] = exp(-(pow(row_idx, 2) + pow(col_idx, 2)) / (2 * pow(sigma, 2)));
            sum += kernel[row_idx + half_size][col_idx + half_size];
        }
    }

    //normalizing
    for (int row_idx = 0; row_idx < size; ++row_idx) {
        for (int col_idx = 0; col_idx < size; ++col_idx) {
            kernel[row_idx][col_idx] /= sum;
        }
    }
    return kernel;
}

QImage apply_filter(const QImage& input_image, const std::vector<std::vector<double>>& kernel) {
    int width = input_image.width();
    int height = input_image.height();
    int half_size = kernel.size() / 2;
    QImage output_image(width, height, QImage::Format_Grayscale8);

    for (int y = half_size; y < height - half_size; ++y) {
        for (int x = half_size; x < width - half_size; ++x) {
            double sum = 0.0;

            for (int kern_row_idx = -half_size; kern_row_idx <= half_size; ++kern_row_idx) {
                for (int kern_col_idx = -half_size; kern_col_idx <= half_size; ++kern_col_idx) {
                    sum += kernel[kern_row_idx + half_size][kern_col_idx + half_size] *
                           //since the RGB values are the same, we get one of them to get the intensity of the grayscale pixel
                           qRed(input_image.pixel(x + kern_col_idx, y + kern_row_idx));
                }
            }
            sum = std::clamp(sum, INTENSITY_MIN, INTENSITY_MAX); //bound x within 0 and 255
            output_image.setPixel(x, y, qRgb(sum, sum, sum)); //a grayscale pixel has RGB values of equal intensity
        }
    }
    return output_image;
}
