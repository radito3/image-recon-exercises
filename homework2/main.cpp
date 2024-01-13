#include <QApplication>
#include <QImage>
#include <QRgb>
#include <QLabel>
#include <QTextStream>
#include <QtGui>
#include <QString>

#include <algorithm>
#include <vector>
#include <iostream>

QImage applyAdaptiveMedianFilter(const QImage& inputImage, int maxWindowSize);
int getWindowSize(const QImage& image, int x, int y, int maxWindowSize);
int getMedian(const QImage& image, int x, int y, int windowSize);
int getPixel(const QImage& image, int x, int y);

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "no file name provided" << std::endl;
        std::cerr << "Usage: <program> <path/to/grayscale_image>" << std::endl;
        return 1;
    }

    QString fileName(QString::fromLocal8Bit(argv[1]));
    QApplication app(argc, argv);
    QImage image;
    QLabel label;

    if (!image.load(fileName)) {
        QTextStream(stderr) << "Failed to load the image." << Qt::endl;
        return -1;
    }

    int maxWindowSize = 7;
    QImage result = applyAdaptiveMedianFilter(image, maxWindowSize);

    label.setPixmap(QPixmap::fromImage(result));
    label.show();

    return app.exec();
}

QImage applyAdaptiveMedianFilter(const QImage& inputImage, int maxWindowSize) {
    QImage outputImage(inputImage);

    for (int y = 0; y < inputImage.height(); ++y) {
        for (int x = 0; x < inputImage.width(); ++x) {
            int windowSize = getWindowSize(inputImage, x, y, maxWindowSize);
            int median = getMedian(inputImage, x, y, windowSize);
            outputImage.setPixel(x, y, qRgb(median, median, median));
        }
    }

    return outputImage;
}

int getWindowSize(const QImage& image, int x, int y, int maxWindowSize) {
    int windowSize = 3; //initial

    while (windowSize <= maxWindowSize) {
        std::vector<int> values;

        for (int i = -windowSize / 2; i <= windowSize / 2; ++i) {
            for (int j = -windowSize / 2; j <= windowSize / 2; ++j) {
                int newX = std::clamp(x + i, 0, image.width() - 1);
                int newY = std::clamp(y + j, 0, image.height() - 1);
                values.push_back(getPixel(image, newX, newY));
            }
        }

        std::sort(values.begin(), values.end());
        int minValue = values.front();
        int maxValue = values.back();
        int medianValue = values[values.size() / 2];

        if (minValue < medianValue && medianValue < maxValue) {
            return windowSize;
        }
        ++windowSize;
    }

    return maxWindowSize;
}

int getMedian(const QImage& image, int x, int y, int windowSize) {
    std::vector<int> values;

    for (int i = -windowSize / 2; i <= windowSize / 2; ++i) {
        for (int j = -windowSize / 2; j <= windowSize / 2; ++j) {
            int newX = std::clamp(x + i, 0, image.width() - 1);
            int newY = std::clamp(y + j, 0, image.height() - 1);
            values.push_back(getPixel(image, newX, newY));
        }
    }

    std::sort(values.begin(), values.end());
    int medianValue = values[values.size() / 2];

    // int minValue = values.front();
    // int maxValue = values.back();
    // int pixelAtXY = getPixel(image, x, y);

    // if (minValue < pixelAtXY && pixelAtXY < maxValue) {
        // this produces worse results at image reconstruction than only returning the median
    //     return pixelAtXY;
    // }
    return medianValue;
}

int getPixel(const QImage& image, int x, int y) {
    return qGray(image.pixel(x, y));
}
