#include <QString>
#include <QtGui>
#include <QApplication>
#include <QImage>
#include <QLabel>
#include <QTextStream>

#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>

QImage apply_filter(const QImage&, int);

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "file name and filter size not provided" << std::endl;
        std::cerr << "Usage: <program> <path/to/grayscale_image> <filter_size>" << std::endl;
        return 1;
    }

    QString file_name(QString::fromLocal8Bit(argv[1]));
    auto filter_size = strtoul(argv[2], NULL, 10);
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

    QImage out_image = apply_filter(image, static_cast<int>(filter_size));

    label.setPixmap(QPixmap::fromImage(out_image));
    label.show();

    return app.exec();
}

QImage apply_filter(const QImage& input_image, int filter_size) {
    int width = input_image.width();
    int height = input_image.height();
    int half_filter_size = filter_size / 2;

    QImage output_image(width, height, QImage::Format_Grayscale8);

    for (int y = 0; y < height; ++y) {
        int start_x = 0;
        int end_x = width - 1;
        int step_x = 1;

        if (y % 2 != 0) {
            start_x = width - 1;
            end_x = 0;
            step_x = -1;
        }

        for (int x = start_x; x != end_x + step_x; x += step_x) {
            std::vector<int> values;
            values.reserve(filter_size);
    
            for (int i = -half_filter_size; i <= half_filter_size; ++i) {
                //set the lower and upper bounds for column value px
                int px = std::clamp(x + i, 0, width - 1);
                //since int is very small and cheap to create, it's irrelevant whether we use push_back() or emplace_back()
                values.push_back(qRed(input_image.pixel(px, y)));
            }

            std::sort(values.begin(), values.end());

            int median;
            if (filter_size % 2 == 0) {
                //vector.at() has bounds check, where as the [] operator does not
                median = (values.at(filter_size / 2 - 1) + values.at(filter_size / 2)) / 2;
            } else {
                median = values.at(filter_size / 2);
            }
            //since the image is grayscale, all RGB values are the same - the intensity of the pixel
            output_image.setPixel(x, y, qRgb(median, median, median));
        }
    }
    return output_image;
}
