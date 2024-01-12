#include <QString>
#include <QtGui>
#include <QApplication>
#include <QImage>
#include <QLabel>
#include <QPainter>
#include <QTextStream>

#include <cmath>
#include <algorithm>
#include <iostream>

const QString file_name = "gray_lena.jpeg";

double contraharmonicMean(const QImage &image, int x, int y, int windowSize, double order) {
    double numerator = 0.0;
    double denominator = 0.0;

    int width = image.width();
    int height = image.height();
  
    int halfWindowSize = windowSize / 2;

    for (int i = -halfWindowSize; i <= halfWindowSize; ++i) {
        for (int j = -halfWindowSize; j <= halfWindowSize; ++j) {
            int posX = std::clamp(x + i, 0, width - 1);
            int posY = std::clamp(y + j, 0, height - 1);

            uchar intensity = image.pixelColor(posX, posY).red();

            numerator += std::pow(intensity, order + 1);
            denominator += std::pow(intensity, order);
        }
    }

    return numerator / denominator;
}

QImage applyContraharmonicMeanFilter(const QImage& image, int windowSize, double order) {
    QImage resultImage(image);

    for (int x = 0; x < image.width(); ++x) {
        for (int y = 0; y < image.height(); ++y) {
            double filteredValue = contraharmonicMean(image, x, y, windowSize, order);
            resultImage.setPixelColor(x, y, qRgb(filteredValue, filteredValue, filteredValue));
        }
    }

    return resultImage;
}

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);
    QImage image;
    QLabel label;

    if (!image.load(file_name)) {
        QTextStream(stderr) << "Failed to load image: " << file_name << Qt::endl;
        return 1;
    }

    int filterSize = 3;
    double Q; // order
    std::cin >> Q;

    QImage result = applyContraharmonicMeanFilter(image, filterSize, Q);

    label.setPixmap(QPixmap::fromImage(result));
    label.show();

    return app.exec();
}
