#include <QString>
#include <QtGui>
#include <QApplication>
#include <QImage>
#include <QLabel>
#include <QPainter>
#include <QTextStream>

#include <string>
#include <vector>
#include <complex>
#include <cmath>

// using namespace std;

const QString file_name = "blurred-lena.jpg";

const int KERN_DIM = 3;

void to_gray(QImage& image) {
    //tranformation coefficients from color to grayscale
    const double red_coef = 0.2989,
                 green_coef = 0.5870,
                 blue_coef = 0.1140;
    
    for (int row_i = 0; row_i < image.height(); row_i++) {
        QRgb* row_p = (QRgb*) image.scanLine(row_i);

        for (int col_i = 0; col_i < image.width(); col_i++) {
            auto pixel = row_p[col_i];
            int red = qRed(pixel);
            int green = qGreen(pixel);
            int blue = qBlue(pixel);

            int gray = red_coef * red + green_coef * green + blue_coef * blue;
            row_p[col_i] = qRgb(gray, gray, gray);
        }
    }
}

QImage arithmetic_mean_filter(const QImage& image, int kern_dim = KERN_DIM) {
    QImage result(image.width(), image.height(), QImage::Format_Grayscale8);
    int half_filter_size = kern_dim / 2;
    double coef = 1.0 / (kern_dim * kern_dim);

    for (int row_idx = 0; row_idx < image.height(); row_idx++) {
        quint8* out_row_ptr = (quint8*) (result.bits() + row_idx * result.bytesPerLine());
        
        for (int col_idx = 0; col_idx < result.width(); col_idx++) {
            double conv_value = 0;
            
            for (int ker_row_idx = 0; ker_row_idx < kern_dim; ker_row_idx++) {
                int x = row_idx - half_filter_size + ker_row_idx;
            
                if (x >= 0 && x < image.height()) {
                    quint8* in_row_ptr = (quint8*) (image.bits() + x * image.bytesPerLine());
                    
                    for (int ker_col_idx = 0; ker_col_idx < kern_dim; ker_col_idx++) {
                        int y = col_idx - half_filter_size + ker_col_idx;
                        
                        if (y >= 0 && y < image.width()) {
                            conv_value += in_row_ptr[y] * coef;
                        }
                    }
                }
            }
            out_row_ptr[col_idx] = (int) conv_value;
        }
    }

    return result;
}

QImage geometric_mean_filter(const QImage& image, int kern_dim = KERN_DIM) {
    QImage result(image.width(), image.height(), QImage::Format_Grayscale8);
    int half_filter_size = kern_dim / 2;
    double coef = 1.0 / (kern_dim * kern_dim);

    for (int row_idx = 0; row_idx < image.height(); row_idx++) {
        quint8* out_row_ptr = (quint8*) (result.bits() + row_idx * result.bytesPerLine());
        
        for (int col_idx = 0; col_idx < result.width(); col_idx++) {
            double prod = 1;
            
            for (int ker_row_idx = 0; ker_row_idx < kern_dim; ker_row_idx++) {
                int x = row_idx - half_filter_size + ker_row_idx;
            
                if (x >= 0 && x < image.height()) {
                    quint8* in_row_ptr = (quint8*) (image.bits() + x * image.bytesPerLine());
                    
                    for (int ker_col_idx = 0; ker_col_idx < kern_dim; ker_col_idx++) {
                        int y = col_idx - half_filter_size + ker_col_idx;
                        
                        if (y >= 0 && y < image.width()) {
                            prod *= in_row_ptr[y];
                        }
                    }
                }
            }
            out_row_ptr[col_idx] = (int) std::pow(prod, coef);
        }
    }

    return result;
}

QImage harmonic_mean_filter(const QImage& image, int kern_dim = KERN_DIM) {
    QImage result(image.width(), image.height(), QImage::Format_Grayscale8);
    int half_filter_size = kern_dim / 2;

    for (int row_idx = 0; row_idx < image.height(); row_idx++) {
        quint8* out_row_ptr = (quint8*) (result.bits() + row_idx * result.bytesPerLine());
        
        for (int col_idx = 0; col_idx < result.width(); col_idx++) {
            double sum = 0;
            
            for (int ker_row_idx = 0; ker_row_idx < kern_dim; ker_row_idx++) {
                int x = row_idx - half_filter_size + ker_row_idx;
            
                if (x >= 0 && x < image.height()) {
                    quint8* in_row_ptr = (quint8*) (image.bits() + x * image.bytesPerLine());
                    
                    for (int ker_col_idx = 0; ker_col_idx < kern_dim; ker_col_idx++) {
                        int y = col_idx - half_filter_size + ker_col_idx;
                        
                        if (y >= 0 && y < image.width()) {
                            sum += 1.0 / in_row_ptr[y];
                        }
                    }
                }
            }
            out_row_ptr[col_idx] = (int) ((kern_dim * kern_dim) / sum);
        }
    }

    return result;
}

int median(std::vector<int>& intens) {
    std::sort(intens.begin(), intens.end());
    intens.erase(std::unique(intens.begin(), intens.end()), intens.end());
    return intens[intens.size() / 2];
}

//salt and pepper filter
QImage median_filter(const QImage& image, int kern_dim = KERN_DIM) {
    QImage result(image.width(), image.height(), QImage::Format_Grayscale8);
    int half_filter_size = kern_dim / 2;

    for (int row_idx = 0; row_idx < image.height(); row_idx++) {
        quint8* out_row_ptr = (quint8*) (result.bits() + row_idx * result.bytesPerLine());
        
        for (int col_idx = 0; col_idx < result.width(); col_idx++) {
            std::vector<int> intens;
            
            for (int ker_row_idx = 0; ker_row_idx < kern_dim; ker_row_idx++) {
                int x = row_idx - half_filter_size + ker_row_idx;
            
                if (x >= 0 && x < image.height()) {
                    quint8* in_row_ptr = (quint8*) (image.bits() + x * image.bytesPerLine());
                    
                    for (int ker_col_idx = 0; ker_col_idx < kern_dim; ker_col_idx++) {
                        int y = col_idx - half_filter_size + ker_col_idx;
                        
                        if (y >= 0 && y < image.width()) {
                            intens.push_back(in_row_ptr[y]);
                        }
                    }
                }
            }
            out_row_ptr[col_idx] = median(intens);
        }
    }

    return result;
}

int main(int argc, char* argv[]) {
    QApplication app(argc, argv);
    QImage image;
    QLabel label;

    if (argc < 2) {
        QTextStream(stderr) << "not enough arguments. need 1" << Qt::endl;
        return 1;
    }

    if (!image.load(file_name)) {
        QTextStream(stderr) << "Failed to load image: " << file_name << Qt::endl;
        return 1;
    }

    to_gray(image);
    image = image.convertToFormat(QImage::Format_Grayscale8);

    QImage out_image;
    if (std::string(argv[1]) == "amean") {
        out_image = arithmetic_mean_filter(image);
    } else if (std::string(argv[1]) == "gmean") {
        out_image = geometric_mean_filter(image);
    } else if (std::string(argv[1]) == "hmean") {
        out_image = harmonic_mean_filter(image);
    } else if (std::string(argv[1]) == "median") {
        //the higher the kernel dimension, the more the noise is removed but also the more blurry the image gets 
        out_image = median_filter(image, 5);
    } else {
        QTextStream(stderr) << "bruh" << Qt::endl;
        return 1;
    }

    label.setPixmap(QPixmap::fromImage(out_image));
    label.show();

    return app.exec();
}
