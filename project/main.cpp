#include <QApplication>
#include <QImage>
#include <QColor>
#include <QRandomGenerator>
#include <QMouseEvent>
#include <QPixmap>
#include <QPainter>
#include <QFileDialog>
#include <QDebug>
#include <QLabel>
#include <QVBoxLayout>
#include <QPushButton>
#include <QWidget>
#include <QRubberBand>

#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <complex>

// ----------------- edge detection -----------------

// Function to apply Gaussian blur to the image
QImage applyGaussianBlur(const QImage &inputImage, int kernelSize, double sigma) {
    QImage outputImage = inputImage;

    int halfKernelSize = kernelSize / 2;
    int width = inputImage.width();
    int height = inputImage.height();

    // Create a 2D kernel for Gaussian blur
    double **kernel = new double*[kernelSize];
    for (int i = 0; i < kernelSize; ++i) {
        kernel[i] = new double[kernelSize];
    }

    // Populate the kernel values based on Gaussian function
    for (int i = -halfKernelSize; i <= halfKernelSize; ++i) {
        for (int j = -halfKernelSize; j <= halfKernelSize; ++j) {
            kernel[i + halfKernelSize][j + halfKernelSize] =
                exp(-(i * i + j * j) / (2 * sigma * sigma)) / (2 * M_PI * sigma * sigma);
        }
    }

    // Normalize the kernel
    double sum = 0.0;
    for (int i = 0; i < kernelSize; ++i) {
        for (int j = 0; j < kernelSize; ++j) {
            sum += kernel[i][j];
        }
    }

    for (int i = 0; i < kernelSize; ++i) {
        for (int j = 0; j < kernelSize; ++j) {
            kernel[i][j] /= sum;
        }
    }

    // Apply the Gaussian blur to the image
    for (int y = halfKernelSize; y < height - halfKernelSize; ++y) {
        for (int x = halfKernelSize; x < width - halfKernelSize; ++x) {
            double sumR = 0.0, sumG = 0.0, sumB = 0.0;

            for (int i = -halfKernelSize; i <= halfKernelSize; ++i) {
                for (int j = -halfKernelSize; i <= halfKernelSize; ++i) {
                    QRgb pixel = inputImage.pixel(x + j, y + i);
                    sumR += qRed(pixel) * kernel[i + halfKernelSize][j + halfKernelSize];
                    sumG += qGreen(pixel) * kernel[i + halfKernelSize][j + halfKernelSize];
                    sumB += qBlue(pixel) * kernel[i + halfKernelSize][j + halfKernelSize];
                }
            }

            outputImage.setPixel(x, y, qRgb(static_cast<int>(sumR), static_cast<int>(sumG), static_cast<int>(sumB)));
        }
    }

    // Clean up memory
    for (int i = 0; i < kernelSize; ++i) {
        delete[] kernel[i];
    }
    delete[] kernel;

    return outputImage;
}

// Function to compute gradients using Sobel operators
void computeGradients(const QImage &inputImage, QImage &gradientMagnitude, QImage &gradientDirection) {
    int width = inputImage.width();
    int height = inputImage.height();

    gradientMagnitude = QImage(width, height, QImage::Format_Grayscale8);
    gradientDirection = QImage(width, height, QImage::Format_Grayscale8);

    for (int y = 1; y < height - 1; ++y) {
        for (int x = 1; x < width - 1; ++x) {
            int gx = qGray(inputImage.pixel(x + 1, y - 1)) + 2 * qGray(inputImage.pixel(x + 1, y)) + qGray(inputImage.pixel(x + 1, y + 1))
                     - qGray(inputImage.pixel(x - 1, y - 1)) - 2 * qGray(inputImage.pixel(x - 1, y)) - qGray(inputImage.pixel(x - 1, y + 1));

            int gy = qGray(inputImage.pixel(x - 1, y + 1)) + 2 * qGray(inputImage.pixel(x, y + 1)) + qGray(inputImage.pixel(x + 1, y + 1))
                     - qGray(inputImage.pixel(x - 1, y - 1)) - 2 * qGray(inputImage.pixel(x, y - 1)) - qGray(inputImage.pixel(x + 1, y - 1));

            int magnitude = static_cast<int>(sqrt(gx * gx + gy * gy));
            gradientMagnitude.setPixel(x, y, qRgb(magnitude, magnitude, magnitude));

            double direction = atan2(gy, gx) * 180.0 / M_PI;
            gradientDirection.setPixel(x, y, qRgb(static_cast<int>(direction), static_cast<int>(direction), static_cast<int>(direction)));
        }
    }
}

// Function for non-maximum suppression
void nonMaximumSuppression(const QImage &gradientMagnitude, const QImage &gradientDirection, QImage &outputImage) {
    int width = gradientMagnitude.width();
    int height = gradientMagnitude.height();

    outputImage = QImage(width, height, QImage::Format_Grayscale8);

    for (int y = 1; y < height - 1; ++y) {
        for (int x = 1; x < width - 1; ++x) {
            int mag = qRed(gradientMagnitude.pixel(x, y));

            int angle = qRed(gradientDirection.pixel(x, y));
            angle = (angle < 0) ? angle + 180 : angle;

            int pixel1, pixel2;

            // Determine the pixels to compare based on the gradient direction
            if ((angle >= 0 && angle < 22.5) || (angle >= 157.5 && angle <= 180)) {
                pixel1 = qRed(gradientMagnitude.pixel(x + 1, y));
                pixel2 = qRed(gradientMagnitude.pixel(x - 1, y));
            } else if (angle >= 22.5 && angle < 67.5) {
                pixel1 = qRed(gradientMagnitude.pixel(x + 1, y - 1));
                pixel2 = qRed(gradientMagnitude.pixel(x - 1, y + 1));
            } else if (angle >= 67.5 && angle < 112.5) {
                pixel1 = qRed(gradientMagnitude.pixel(x, y - 1));
                pixel2 = qRed(gradientMagnitude.pixel(x, y + 1));
            } else {  // angle >= 112.5 && angle < 157.5
                pixel1 = qRed(gradientMagnitude.pixel(x - 1, y - 1));
                pixel2 = qRed(gradientMagnitude.pixel(x + 1, y + 1));
            }

            // Perform non-maximum suppression
            if (mag >= pixel1 && mag >= pixel2) {
                outputImage.setPixel(x, y, qRgb(mag, mag, mag));
            } else {
                outputImage.setPixel(x, y, qRgb(0, 0, 0));
            }
        }
    }
}

// Function for edge tracking by hysteresis
void edgeTrackingByHysteresis(const QImage &inputImage, QImage &outputImage, double lowThreshold, double highThreshold) {
    int width = inputImage.width();
    int height = inputImage.height();

    outputImage = QImage(width, height, QImage::Format_Grayscale8);

    // Define weak and strong edge intensity values
    int weakEdge = 50;
    int strongEdge = 255;

    // Initialize the output image with weak edges
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            if (qRed(inputImage.pixel(x, y)) > lowThreshold) {
                outputImage.setPixel(x, y, qRgb(weakEdge, weakEdge, weakEdge));
            } else {
                outputImage.setPixel(x, y, qRgb(0, 0, 0));
            }
        }
    }

    // Identify strong edges based on the high threshold
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            if (qRed(inputImage.pixel(x, y)) > highThreshold) {
                outputImage.setPixel(x, y, qRgb(strongEdge, strongEdge, strongEdge));
            }
        }
    }

    // Perform edge tracking by hysteresis
    for (int y = 1; y < height - 1; ++y) {
        for (int x = 1; x < width - 1; ++x) {
            if (qRed(outputImage.pixel(x, y)) == weakEdge) {
                // Check neighboring pixels for strong edges
                for (int i = -1; i <= 1; ++i) {
                    for (int j = -1; j <= 1; ++j) {
                        if (qRed(outputImage.pixel(x + j, y + i)) == strongEdge) {
                            outputImage.setPixel(x, y, qRgb(strongEdge, strongEdge, strongEdge));
                            break;
                        }
                    }
                }
            }
        }
    }
}

// ----------------- edge detection -----------------

// ---- util ----

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

// ---- util ----

// ----------------- binarization -----------------

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
    }
}

// ----------------- binarization -----------------

// ----------------- segmentation -----------------

void depthFirstSearch(const QImage &binaryImage, std::vector<std::vector<int>> &labeledImage, int x, int y, int label) {
    if (x < 0 || x >= binaryImage.width() || y < 0 || y >= binaryImage.height())
        return;

    if (binaryImage.pixelColor(x, y).black() == 1 && labeledImage[y][x] == 0) {
        labeledImage[y][x] = label;

        depthFirstSearch(binaryImage, labeledImage, x - 1, y - 1, label);
        depthFirstSearch(binaryImage, labeledImage, x, y - 1, label);
        depthFirstSearch(binaryImage, labeledImage, x + 1, y - 1, label);
        depthFirstSearch(binaryImage, labeledImage, x - 1, y, label);
        depthFirstSearch(binaryImage, labeledImage, x + 1, y, label);
        depthFirstSearch(binaryImage, labeledImage, x - 1, y + 1, label);
        depthFirstSearch(binaryImage, labeledImage, x, y + 1, label);
        depthFirstSearch(binaryImage, labeledImage, x + 1, y + 1, label);
    }
}

int connectedComponentAnalysis(const QImage &binaryImage, std::vector<std::vector<int>> &labeledImage) {
    labeledImage.resize(binaryImage.height(), std::vector<int>(binaryImage.width(), 0));
    int label = 0;

    for (int y = 0; y < binaryImage.height(); ++y) {
        for (int x = 0; x < binaryImage.width(); ++x) {
            if (binaryImage.pixelColor(x, y).black() == 1 && labeledImage[y][x] == 0) {
                label++;
                depthFirstSearch(binaryImage, labeledImage, x, y, label);
            }
        }
    }

    return label;
}

QImage generateSegmentedImage(const QImage &binaryImage, const std::vector<std::vector<int>> &labeledImage, int numLabels) {
    QImage segmented(binaryImage.size(), QImage::Format_RGB32);
    std::vector<QColor> labelColors(numLabels);

    for (int i = 0; i < numLabels; ++i) {
        labelColors[i] = QColor::fromRgb(QRandomGenerator::global()->bounded(0, 256) % 256,
                                         QRandomGenerator::global()->bounded(0, 256) % 256, 
                                         QRandomGenerator::global()->bounded(0, 256) % 256);
    }

    for (int y = 0; y < binaryImage.height(); ++y) {
        for (int x = 0; x < binaryImage.width(); ++x) {
            int label = labeledImage[y][x];
            QColor color = label > 0 ? labelColors[label - 1] : QColor(Qt::white); // Set background to white
            segmented.setPixelColor(x, y, color);
        }
    }

    return segmented;
}

QImage segmentForegroundBackground(const QImage &inputImage, int threshold = 128) {
    int width = inputImage.width();
    int height = inputImage.height();

    QImage outputImage = QImage(width, height, QImage::Format_Mono);

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            // Get the pixel color (0 or 1 for binary image)
            QRgb pixelColor = inputImage.pixel(x, y);
            int value = qRed(pixelColor);

            // Apply thresholding to segment foreground from background
            if (value > threshold) {
                outputImage.setPixel(x, y, 1); // Foreground
            } else {
                outputImage.setPixel(x, y, 0); // Background
            }
        }
    }
    return outputImage;
}

// ----------------- segmentation -----------------

class ImageViewer : public QWidget {
    Q_OBJECT

public:
    ImageViewer(QWidget *parent = nullptr) 
        : QWidget(parent), isImageDisplayed(false), isSelectingRect(false)
    {
        setupUi();
    }

protected:
    void mousePressEvent(QMouseEvent *event) override {
        if (isSelectingRect && event->button() == Qt::LeftButton) {
            selectionRect.setTopLeft(event->pos());

            origin = event->pos();
            if (!rubberBand)
                rubberBand = new QRubberBand(QRubberBand::Rectangle, this);
            rubberBand->setGeometry(QRect(origin, QSize()));
            rubberBand->show();
        }
    }

    void mouseMoveEvent(QMouseEvent *event) override {
        if (isSelectingRect) {
            selectionRect.setBottomRight(event->pos());

            rubberBand->setGeometry(QRect(origin, event->pos()).normalized());
        }
    }

    void mouseReleaseEvent(QMouseEvent *event) override {
        if (isSelectingRect && event->button() == Qt::LeftButton) {
            rubberBand->hide();

            finishCrop();
        }
    }

private slots:
    void openImage() {
        QString imagePath = QFileDialog::getOpenFileName(this, "Open Image", "", "Images (*.png *.jpg *.bmp)");

        if (!imagePath.isEmpty()) {
            displayImage(imagePath);
        }
    }

    void startCrop() {
        isSelectingRect = true;
        imageLabel->setCursor(Qt::CrossCursor);
    }

    void finishCrop() {
        isSelectingRect = false;
        imageLabel->setCursor(Qt::ArrowCursor);

        // Ensure the region of interest QRect is within the image bounds
        QRect roiRect = selectionRect.intersected(image.rect());

        QImage procssingImg = image.copy(roiRect);

        QImage binaryImg = transformImage(procssingImg);

        // TODO: maybe apply a gaussian filter before contour detection as it is a bit inaccurate
        auto contours = detectContours(binaryImg);

        // draw contours over the original image
        QPainter painter(&image);
        painter.setPen(Qt::red);

        for (const auto& contour : contours) {
            double x = std::real(contour);
            double y = std::imag(contour);

            // Draw points on the original image within the specified QRect
            painter.drawPoint(roiRect.x() + static_cast<int>(x), roiRect.y() + static_cast<int>(y));
        }

        imageLabel->setPixmap(QPixmap::fromImage(image));
    }

private:
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
                        contours.push_back(std::complex<double>(x, y));
                    }
                }
            }
        }
        return contours;
    }

    void setupUi() {
        imageLabel = new QLabel(this);
        imageLabel->setScaledContents(true);

        layout = new QVBoxLayout(this);
        layout->addWidget(imageLabel);

        openButton = new QPushButton("Open Image", this);
        layout->addWidget(openButton);

        connect(openButton, &QPushButton::clicked, this, &ImageViewer::openImage);
    }

    void displayImage(const QString &imagePath) {
        QPixmap image(imagePath);
        if (image.isNull()) {
            qDebug() << "Failed to load image:" << imagePath;
            return;
        }

        this->image = image.toImage();

        imageLabel->setPixmap(image);
        setWindowTitle("Image Viewer - " + imagePath);

        if (!isImageDisplayed) {
            isImageDisplayed = true;
            cropButton = new QPushButton("Select region", this);
            layout->addWidget(cropButton);

            connect(cropButton, &QPushButton::clicked, this, &ImageViewer::startCrop);

            // TODO: add button for segmentation of foreground and background
        }
    }

    QImage transformImage(const QImage& image) {
        QImage copy = image;
        if (!image.allGray()) {
            to_gray(copy);
        }
        QImage grayImage = copy.convertToFormat(QImage::Format_Grayscale8);
        toBinImage(grayImage);
        return grayImage;
        // return segmentForegroundBackground(grayImage);

        // std::vector<std::vector<int>> labeledImage;
        // int numLabels = connectedComponentAnalysis(grayImage, labeledImage);
        // return generateSegmentedImage(grayImage, labeledImage, numLabels);
    }

    QPoint origin;
    QRubberBand* rubberBand;

    QImage image;
    QLabel *imageLabel;
    QVBoxLayout *layout;
    QPushButton *openButton;

    bool isImageDisplayed;
    QPushButton *cropButton;
    bool isSelectingRect;
    QRect selectionRect;
};

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);

    // first 3 steps are uneccessary for the MNIST dataset, as the images are very small and do not have much noise

    // Step 1: Apply Gaussian blur
    // QImage blurredImage = applyGaussianBlur(inputImage, 5, 1.0);

    // Step 2: Compute gradients
    // QImage gradientMagnitude, gradientDirection;
    // computeGradients(inputImage, gradientMagnitude, gradientDirection);

    // Step 3: Non-maximum suppression
    // QImage nonMaxSuppression;
    // nonMaximumSuppression(gradientMagnitude, gradientDirection, nonMaxSuppression);

    // Step 4: Edge tracking by hysteresis
    // QImage finalEdges;
    // edgeTrackingByHysteresis(inputImage, /*nonMaxSuppression,*/ finalEdges, 70.0, 200.0);

    // finalEdges.save("output.jpg");

    ImageViewer viewer;
    viewer.show();

    return app.exec();
}

#include "main.moc"
