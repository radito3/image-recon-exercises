#include <QApplication>
#include <QWidget>
#include <QLabel>
#include <QFile>
#include <QDebug>
#include <QImage>

#include <iostream>
#include <vector>
#include <complex>

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

// very small constant to compare floating-point numbers
const double EPS = 1.0e-14;

// gray-level intensity minimum and maximum
const int INTENS_MIN = 0;
const int INTENS_MAX = 255;

// accumulate normalized histogram
void calcHisto(const QImage &image, double histo[]) {
    for (int i = 0; i <= INTENS_MAX; i++) {
        histo[i] = 0;
    }

    for (int indx_row = 0; indx_row < image.height(); indx_row++) {
        quint8* ptr_row = (quint8*)(image.bits() + indx_row * image.bytesPerLine());
        for (int indx_col = 0; indx_col < image.width(); indx_col++) {
            histo[ptr_row[indx_col]]++;
        }
    }

    int numb_pix = image.height() * image.width();
    for (int i = 0; i <= INTENS_MAX; i++) {
        histo[i] /= numb_pix;
    }
}

// apply threshold on a gray-scale image to calculate binary image
void thresh(QImage &image, int thr) {
    for (int indx_row = 0; indx_row < image.height(); indx_row++) {
        quint8* ptr_row = (quint8*)(image.bits() + indx_row * image.bytesPerLine());
        for (int indx_col = 0; indx_col < image.width(); indx_col++) {
            ptr_row[indx_col] = (ptr_row[indx_col] < thr) ? INTENS_MIN : INTENS_MAX;
        }
    }
}

// Otsu's method
int otsu(const double histo[]) {
    // compute cumulative sums
    double p_1[INTENS_MAX + 1] = {0};
    p_1[0] = histo[0];
    for (int i = 1; i <= INTENS_MAX; i++)
    {
        p_1[i] = p_1[i - 1] + histo[i];
    }

    // cumulative mean
    double m[INTENS_MAX + 1] = {0};
    for (int i = 1; i <= INTENS_MAX; i++)
    {
        m[i] = m[i - 1] + i * histo[i];
    }

    // global mean
    double m_g = m[INTENS_MAX];

    // between-class
    double b_c[INTENS_MAX + 1] = {0};
    for (int i = 1; i <= INTENS_MAX; i++)
    {
        double div = (p_1[i] * (1 - p_1[i]));
        b_c[i] = 
            fabs(div < EPS) ? 0 :
            ((m_g * p_1[i] - m[i]) * (m_g * p_1[i] - m[i])) / div;
    }

    // find max
    double max = 0;
    int max_i = 0;
    for (int i = 0; i <= INTENS_MAX; i++)
    {
        if (b_c[i] > max)
        {
            max = b_c[i];
            max_i = i;
        }
    }

    return max_i;
}

// convert a gray-scale image to binary image
void toBinImage(QImage &image) {
    if (image.format() == QImage::Format_Grayscale8) {
        double histo[INTENS_MAX + 1];
        calcHisto(image, histo);  
        int th = otsu(histo);
        thresh(image, th);
    } else {
        std::cout << "not grayscale" << std::endl;
    }
}

std::vector<std::complex<double>> detectContours(const QImage& binaryImage) {
    std::vector<std::complex<double>> contours;

    for (int y = 0; y < binaryImage.height(); ++y) {
        for (int x = 0; x < binaryImage.width(); ++x) {
            // Check if the pixel is part of the object (black pixel in binary image)
            if (qRed(binaryImage.pixel(x, y)) == 0) {
                // Check neighbors to determine if it is on the boundary
                bool isBoundary = false;

                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dx = -1; dx <= 1; ++dx) {
                        int nx = x + dx;
                        int ny = y + dy;

                        // Skip pixels outside the image boundary
                        if (nx >= 0 && nx < binaryImage.width() && ny >= 0 && ny < binaryImage.height()) {
                            // Check if the neighbor is outside the object (white pixel in binary image)
                            if (qRed(binaryImage.pixel(nx, ny)) == 255) {
                                isBoundary = true;
                                break;
                            }
                        }
                    }

                    if (isBoundary) break;
                }

                // If the pixel is on the boundary, add it to the contours
                if (isBoundary) {
                    contours.emplace_back(x, y);
                }
            }
        }
    }
    return contours;
}

void fft(std::vector<std::complex<double>>& a) {
    int n = a.size();
    if (n <= 1) {
        return;
    }

    // Divide the vector into even and odd parts
    std::vector<std::complex<double>> even, odd;
    for (int i = 0; i < n; i += 2) {
        even.push_back(a[i]);
        odd.push_back(a[i + 1]);
    }

    // Recursively compute FFT for even and odd parts
    fft(even);
    fft(odd);

    // Combine the results
    for (int i = 0; i < n / 2; i++) {
        // std::polar(r, theta) r - magnitude, theta - phase angle
        // essentially -> std::complex(r * cos(theta), r * sin(theta))
        std::complex<double> t = std::polar(1.0, -2.0 * M_PI * i / n) * odd[i];
        a[i] = even[i] + t;
        a[i + n / 2] = even[i] - t;
    }
}

int main(int argc, char* argv[]) {
    QApplication app(argc, argv);
    QImage image;
    QLabel label;

    if (argc < 2) {
        std::cerr << "not enough arguments. requred 2" << std::endl;
        return 1;
    }

    if (!image.load(QString::fromLocal8Bit(argv[1]))) {
        std::cerr << "could not load image: " << argv[1] << std::endl;
        return 1;
    }

    if (!image.allGray()) {
        to_gray(image);
    }
    image = image.convertToFormat(QImage::Format_Grayscale8);

    // segmenting into object and background
    toBinImage(image);

    std::vector<std::complex<double>> contours = detectContours(image);

    fft(contours);

    for (const auto& transformed : contours) {
        std::cout << transformed;
    }
    std::cout << std::endl;

    // image.save("test_img.png");

    return app.exec();
}
