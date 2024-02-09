/*
* Add the following to the .pro Qt file:
*   QT += gui
*   QT += widgets
*   QT += core network
*   LIBS += -lz
*   
* Optionally:
*   MOC_DIR = $$OUT_PWD/moc
*   OBJECTS_DIR = $$OUT_PWD/obj
*   INCLUDEPATH += $$MOC_DIR
*/
#include <QApplication>
#include <QNetworkAccessManager>
#include <QNetworkRequest>
#include <QNetworkReply>
#include <QProgressBar>
#include <QWidget>
#include <QLabel>
#include <QVBoxLayout>
#include <QProgressDialog>
#include <QFile>
#include <QEventLoop>
#include <QDebug>
// #include <QThreadPool>
#include <QImage>
#include <QPainter>
#include <QPainterPath>
#include <QGraphicsView>
#include <QGraphicsScene>

#include <initializer_list>
#include <unordered_map>
// #include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <bitset>
#include <unordered_map>
#include <complex>
#include <cmath>
#include <cstdio>
#include <limits>
#include <cstdlib>

#include <zlib.h>

class FileDownloader : public QObject {
    Q_OBJECT

public:
    FileDownloader(QUrl&& fileUrl, QProgressBar* progressBar = new QProgressBar, QObject* parent = nullptr) 
        : QObject(parent), url(fileUrl), progressBar(progressBar)
    {
        file.setFileName(url.fileName());
        if (!file.open(QIODevice::WriteOnly)) {
            qDebug() << "Error opening file for writing:" << file.errorString();
        }
        progressBar->setMinimum(0);
        progressBar->setMaximum(0);
    }

    void download() {
        if (!file.isOpen()) {
            qDebug() << "File not open. Can not download...";
            return;
        }

        auto manager = new QNetworkAccessManager(this);

        reply = manager->get(QNetworkRequest(url));

        // Connect signals/slots to handle the data as it arrives
        connect(reply, &QNetworkReply::readyRead, this, &FileDownloader::onReadyRead);
        connect(reply, &QNetworkReply::finished, this, &FileDownloader::onFinished);
        connect(reply, &QNetworkReply::downloadProgress, this, &FileDownloader::onDownloadProgress);

        // Create an event loop to wait for the request to finish
        loop.exec();
    }

    std::filesystem::path getFilePath() const {
        return file.filesystemFileName();
    }

private slots:
    void onReadyRead() {
        while (reply->bytesAvailable()) {
            // Read and write data in chunks of 32 KB
            QByteArray data = reply->read(qMin(reply->bytesAvailable(), qint64(32 * 1024)));
            qint64 bytesWritten = file.write(data);

            // Check for errors
            if (bytesWritten == -1) {
                qDebug() << "Error writing to file:" << file.errorString();
                // You might want to abort the download or handle the error accordingly
                return;
            }
        }
    }

    void onFinished() {
        // Check for errors
        if (reply->error() != QNetworkReply::NoError) {
            qDebug() << "Error downloading" << url.fileName() << ":" << reply->errorString();
        } else {
            qDebug() << "Finished downloading" << url.fileName();
        }

        // Clean up
        file.close();
        reply->deleteLater();
        loop.quit();
    }

    void onDownloadProgress(qint64 bytesReceived, qint64 bytesTotal) {
        if (!progressBar->maximum()) {
            progressBar->setMaximum(static_cast<int>(bytesTotal));
        }
        progressBar->setValue(static_cast<int>(bytesReceived));
    }

private:
    QUrl url;
    QNetworkReply* reply;
    QFile file;
    QProgressBar* progressBar;
    QEventLoop loop;
};

class DownloadDialog : public QWidget {
    Q_OBJECT

public:
    DownloadDialog(std::initializer_list<const char*> urls, QWidget* parent = nullptr) : QWidget(parent) {
        setFixedSize(400, 500);
        setWindowTitle("MNIST Downloader");

        // Create layout
        QVBoxLayout* layout = new QVBoxLayout(this);

        for (auto* url : urls) {
            // Create widgets
            auto pb = new QProgressBar(this);
            pb->setAlignment(Qt::AlignHCenter);

            QLabel* titleLabel = new QLabel(QUrl(url).fileName(), this);
            titleLabel->setAlignment(Qt::AlignHCenter);

            auto pd = new QProgressDialog;
            pd->setWindowModality(Qt::WindowModal);
            pd->setAutoClose(false);
            pd->setAutoReset(false);
            pd->setLabel(titleLabel);
            pd->setBar(pb);
            pd->setCancelButton(nullptr);

            layout->addWidget(pd);

            progressBars.emplace(url, pb);
        }
    }

    std::vector<std::filesystem::path> downloadAll() {
        std::vector<std::filesystem::path> filePaths;
        // downladed sequentially instead of concurrently because the widget
        // does not show up until after all the downloads have finished which
        // invalidates the usage of progress bars in the first place...
        for (auto [url, pb] : progressBars) {
            FileDownloader downloader(QUrl(url), pb);
            // QThreadPool::globalInstance().start([=]() { download(url, pb); });
            downloader.download();
            filePaths.push_back(downloader.getFilePath());
        }
        return filePaths;
    }

private:
    std::unordered_map<const char*, QProgressBar*> progressBars;
};

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

// Function to apply FFT on 1D sequences and return magnitude spectrum
std::vector<double> applyFFT(const std::vector<std::complex<double>> &contours) {
    int N = contours.size();
    std::vector<std::complex<double>> input = contours;

    fft(input);

    // Extract magnitude spectrum
    std::vector<double> spectrum(N);
    for (int i = 0; i < N; i++) {
        // |z| = sqrt(real(z)^2 + imag(z)^2)
        spectrum[i] = std::abs(input[i]);
    }
    return spectrum;
}

// read gzipped MNIST images
std::vector<std::vector<uint8_t>> readMNISTImages(const char* filePath) {
    gzFile file = gzopen(filePath, "rb");

    if (!file) {
        qDebug() << "Error opening file: " << filePath;
        return {};
    }

    // Read magic number and header information
    uint32_t magicNumber, numImages, numRows, numCols;
    gzread(file, &magicNumber, sizeof(magicNumber));
    gzread(file, &numImages, sizeof(numImages));
    gzread(file, &numRows, sizeof(numRows));
    gzread(file, &numCols, sizeof(numCols));

    magicNumber = __builtin_bswap32(magicNumber);
    numImages = __builtin_bswap32(numImages);
    numRows = __builtin_bswap32(numRows);
    numCols = __builtin_bswap32(numCols);

    std::vector<std::vector<uint8_t>> images(numImages, std::vector<uint8_t>(numRows * numCols));

    // Read image data
    for (uint32_t i = 0; i < numImages; ++i) {
        gzread(file, images[i].data(), numRows * numCols);
    }

    gzclose(file);

    return images;
}

// read gzipped MNIST labels
std::vector<uint8_t> readMNISTLabels(const char* filePath) {
    gzFile file = gzopen(filePath, "rb");

    if (!file) {
        qDebug() << "Error opening file:" << filePath;
        return {};
    }

    // Read magic number and header information
    uint32_t magicNumber, numLabels;
    gzread(file, &magicNumber, sizeof(magicNumber));
    gzread(file, &numLabels, sizeof(numLabels));

    magicNumber = __builtin_bswap32(magicNumber);
    numLabels = __builtin_bswap32(numLabels);

    std::vector<uint8_t> labels(numLabels);

    // Read label data
    gzread(file, labels.data(), numLabels);

    gzclose(file);

    return labels;
}

// -------------- contour detection --------------

QPainterPath detectContours(const QImage& binaryImage) {
    QPainterPath contours;

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
                    contours.moveTo(x, y);
                    contours.lineTo(x + 1, y + 1);
                }
            }
        }
    }

    return contours;
}

std::vector<std::complex<double>> detectContoursV(const QImage& binaryImage) {
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

// Alias for a complex number representing a point in 2D space
using Point = std::complex<double>;

// Function to calculate the Euclidean distance between two points (complex numbers)
double distance(const Point& p1, const Point& p2) {
    // sqrt((p2.x - p1.x)^2 + (p2.y - p1.y)^2);
    return std::abs(p2 - p1);
}

// Function to connect nearby contours
void connectContours(std::vector<Point>& contours, double proximityThreshold = 25.0) {
    std::vector<Point> connectedContours;

    for (size_t i = 0; i < contours.size(); ++i) {
        size_t j = (i + 1) % contours.size(); // Circular index to connect the last point with the first
        // Check if the endpoints of the two contours are close enough
        if (distance(contours[i], contours[j]) < proximityThreshold) {
            // Merge the two contours
            connectedContours.insert(connectedContours.end(), contours.begin() + i, contours.begin() + j + 1);
        }
    }

    contours = connectedContours;
}

// -------------- contour detection --------------

// -------------- analysis --------------------

// this might only need be done for the test dataset?
// as in, after the model has been trained:
//  - check its base accuracy
//  - get the high frequencies and start removing them (tune the peak threshold) from contours spectrum
//  - until a certain accuracy regression threshold is reached 
//      (e.g. if the base model is 80% accurate, a 10% reduction is the cut-off point - get the difference in frequencies)
// the result of the experiment would be the number of high frequencies removed and their threshold (in findPeaks) for the
// configured accuracy regression cut-off point (might be a program parameter)

std::vector<size_t> findPeaks(const std::vector<double>& spectrum, double threshold = 0.5) {
    std::vector<size_t> peaks;
    for (size_t i = 1; i < spectrum.size() - 1; ++i) {
        // the threshold is also a hyperparameter -> tune it for best results
        if (spectrum[i] > spectrum[i - 1] && spectrum[i] > spectrum[i + 1] && spectrum[i] > threshold) {
            peaks.push_back(i);
        }
    }
    return peaks; // elements are indicies in the magnitute spectrum vector, pointing to the coefficients
}

/*
    std::cout << "Dominant Frequencies:\n";
    for (size_t peakIndex : peakIndices) {
        // Frequency calculation: peakIndex / N, where N is the size of the spectrum
        double frequency = static_cast<double>(peakIndex) / spectrum.size();
        std::cout << "Frequency: " << frequency << " Hz, Magnitude: " << spectrum[peakIndex] << "\n";
    }
*/

/*
"high" frequencies would correspond to components in the spectrum with higher magnitudes.
 These high-magnitude components represent the presence of rapid changes or variations in the contour.

general guideline for interpreting frequency components in a spectrum:
 - Low frequencies: Correspond to slowly changing or smooth variations in the signal.
 - High frequencies: Correspond to rapid changes or abrupt transitions in the signal.

When analyzing the spectrum of a contour, high-frequency components could capture details such as sharp corners,
 small details, or irregularities in the contour. By identifying and analyzing these high-frequency components,
  you may gain insights into the detailed structure of the contour.

*/

// -------------- analysis --------------------

// --------------- image processing ---------------

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
    }
}

// --------------- image processing ---------------

// --------------- classification ---------------

// 0 % accuracy
class GaussianNaiveBayes {
public:
    void train(const std::vector<std::vector<double>>& X, const std::vector<uint8_t>& y) {
        // Separate data by class
        separateByClass(X, y);

        // Calculate mean and standard deviation for each class and feature
        calculateClassStatistics();
    }

    uint8_t predict(const std::vector<double>& x) {
        double maxProbability = -INFINITY;
        uint8_t predictedClass = 0;

        for (const auto& [currentClass, classStats] : classStatistics) {
            double probability = calculateProbability(x, classStats.mean, classStats.stdev);

            if (probability > maxProbability) {
                maxProbability = probability;
                predictedClass = currentClass;
            }
        }

        return predictedClass;
    }

private:
    struct ClassStatistics {
        std::vector<double> mean;
        std::vector<double> stdev;
    };

    std::unordered_map<uint8_t, ClassStatistics> classStatistics;

    void separateByClass(const std::vector<std::vector<double>>& X, const std::vector<uint8_t>& y) {
        for (size_t i = 0; i < X.size(); ++i) {
            uint8_t label = y[i];
            if (classStatistics.find(label) == classStatistics.end()) {
                classStatistics[label] = ClassStatistics();
            }

            for (size_t j = 0; j < X[i].size(); ++j) {
                classStatistics[label].mean.push_back(X[i][j]);
            }
        }
    }

    void calculateClassStatistics() {
        for (auto& [currentClass, classStats] : classStatistics) {
            qDebug() << "calculating for class" << currentClass;

            size_t numSamples = classStats.mean.size() / classStatistics.size();

            // Calculate mean
            for (double& mean : classStats.mean) {
                mean /= numSamples;
            }

            // Calculate standard deviation
            for (size_t i = 0; i < numSamples; ++i) {
                for (size_t j = 0; j < classStats.mean.size(); ++j) {
                    double diff = classStatistics[currentClass].mean[j] - classStats.mean[j];
                    classStats.stdev.push_back(diff * diff);
                }
            }

            for (double& stdev : classStats.stdev) {
                stdev = sqrt(stdev / numSamples);
            }
        }
    }

    double calculateProbability(const std::vector<double>& x, const std::vector<double>& mean, const std::vector<double>& stdev) {
        double probability = 1.0;
        for (size_t i = 0; i < x.size(); ++i) {
            double exponent = exp(-((x[i] - mean[i]) * (x[i] - mean[i])) / (2 * stdev[i] * stdev[i]));
            probability *= (1 / (sqrt(2 * M_PI) * stdev[i])) * exponent;
        }
        return probability;
    }
};

// 0 % accuracy
class NaiveBayesClassifier {
private:
    // Number of classes
    uint8_t numClasses;

    // Class-wise counts
    std::vector<double> classCounts;

    // Mean and variance per class and feature
    std::vector<std::vector<std::pair<double, double>>> featureStats;

public:
    NaiveBayesClassifier(uint8_t numClasses) : numClasses(numClasses) {
        classCounts.resize(numClasses, 0.0);
        featureStats.resize(numClasses);
    }

    ~NaiveBayesClassifier() {
        // Clear allocated memory
        for (uint8_t c = 0; c < numClasses; ++c) {
            featureStats[c].clear();
        }
    }

    void train(const std::vector<std::vector<double>>& data, const std::vector<uint8_t>& labels) {
        for (size_t i = 0; i < data.size(); ++i) {
            uint8_t label = labels[i];
            classCounts[label] += 1.0;

            if (featureStats[label].empty()) {
                featureStats[label].resize(data[i].size(), {0.0, 0.0});
            }

            for (size_t j = 0; j < data[i].size(); ++j) {
                double featureValue = data[i][j];

                // Update mean and variance (we'll use variance instead of stdev)
                double& mean = featureStats[label][j].first;
                double& variance = featureStats[label][j].second;

                double delta = featureValue - mean;
                mean += delta / classCounts[label];
                variance += delta * (featureValue - mean);
            }
        }

        // Calculate variance and handle cases where variance is zero
        for (uint8_t c = 0; c < numClasses; ++c) {
            for (size_t i = 0; i < featureStats[c].size(); ++i) {
                double& variance = featureStats[c][i].second;
                variance /= classCounts[c];
                if (variance == 0.0) {
                    // Add a small value to avoid division by zero in logProbability
                    variance = std::numeric_limits<double>::epsilon();
                }
            }
        }
    }

    double calculateLogProbability(double x, double mean, double variance) {
        // Log probability density function of normal distribution
        double exponent = -0.5 * pow((x - mean), 2) / variance;
        return -0.5 * log(2 * M_PI * variance) + exponent;
    }

    uint8_t predict(const std::vector<double>& data) {
        double maxLogProbability = -std::numeric_limits<double>::infinity();
        uint8_t predictedClass = 0;

        for (uint8_t c = 0; c < numClasses; ++c) {
            double logProbability = log(classCounts[c] / data.size());  // Class prior

            for (size_t i = 0; i < data.size(); ++i) {
                if (i < featureStats[c].size()) {
                    double mean = featureStats[c][i].first;
                    double variance = featureStats[c][i].second;
                    logProbability += calculateLogProbability(data[i], mean, variance);
                }
            }

            if (logProbability > maxLogProbability) {
                maxLogProbability = logProbability;
                predictedClass = c;
            }
        }

        return predictedClass;
    }
};

// --------------- classification ---------------

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);

    // ------------------- contours experiment -----------------------
    // // Assuming you have a binary image (black and white)
    // QImage binaryImage("bin_img.jpeg");
    // QPainterPath contoursP = detectContours(binaryImage);
    // // TODO: try with and without connecting the contours to see which model performs better
    // std::vector<Point> contours = detectContoursV(binaryImage);
    // // Set the proximity threshold (adjust as needed)
    // // after a certain value (~20-30), it stops affecting the connections
    // double proximityThreshold = 25.0; // this is, essentially, a hyperparameter
    // // Connect nearby contours
    // connectContours(contours, proximityThreshold);

    // // Create a QPainterPath from the connected contours
    // QPainterPath connectedPath;
    // if (!contours.empty()) {
    //     connectedPath.moveTo(contours.front().real(), contours.front().imag());
    //     for (const auto& point : contours) {
    //         connectedPath.lineTo(point.real(), point.imag());
    //     }
    // }

    // // Create a QGraphicsScene and add the contours to it
    // QGraphicsScene scene;
    // scene.addPath(connectedPath);
    // // scene.addPath(contoursP);
    // // Create a QGraphicsView to display the scene
    // QGraphicsView view(&scene);
    // view.show();
    // ------------------- contours experiment -----------------------


    qDebug() << "Downloading MNIST dataset files...";
    DownloadDialog downloadDialog({
        "http://yann.lecun.com/exdb/mnist/train-images-idx3-ubyte.gz", //train images (60 000)
        "http://yann.lecun.com/exdb/mnist/train-labels-idx1-ubyte.gz", //train labels
        "http://yann.lecun.com/exdb/mnist/t10k-images-idx3-ubyte.gz",  //test images (10 000)
        "http://yann.lecun.com/exdb/mnist/t10k-labels-idx1-ubyte.gz"   //test labels
    });
    downloadDialog.show();
    auto filePaths = downloadDialog.downloadAll();

    const char* trainImagesFile = "train-images-idx3-ubyte.gz";
    const char* trainLabelsFile = "train-labels-idx1-ubyte.gz";
    const char* testImagesFile = "t10k-images-idx3-ubyte.gz";
    const char* testLabelsFile = "t10k-labels-idx1-ubyte.gz";

    qDebug() << "Reading MNIST train dataset...";
    std::vector<std::vector<uint8_t>> trainImages = readMNISTImages(trainImagesFile);
    std::vector<uint8_t> trainLabels = readMNISTLabels(trainLabelsFile);
    std::vector<std::vector<uint8_t>> testImages = readMNISTImages(testImagesFile);
    std::vector<uint8_t> testLabels = readMNISTLabels(testLabelsFile);

    qDebug() << "computing image contours and their magnitude spectrums...";
    std::vector<std::vector<double>> trainData(trainImages.size() / 12);
    for (size_t i = 0; i < trainImages.size() / 12; ++i) {
        // Create a QImage from the vector<uint8_t>
        QImage image(trainImages[i].data(), 28, 28, QImage::Format_Grayscale8);
        if (!image.allGray()) {
            to_gray(image);
        }
        QImage grayImage = image.convertToFormat(QImage::Format_Grayscale8);
        toBinImage(grayImage);
        auto contours = detectContoursV(grayImage);
        std::vector<double> magnituteSpectrum = applyFFT(contours);
        trainData[i] = magnituteSpectrum;
    }
  
    qDebug() << "creating classifier";
    GaussianNaiveBayes classifier;

    qDebug() << "training...";
    classifier.train(trainData, trainLabels);

    qDebug() << "testing...";
    std::bitset<10000> guesses;
    for (size_t i = 0; i < testImages.size(); ++i) {
        // Create a QImage from the vector<uint8_t>
        QImage image(testImages[i].data(), 28, 28, QImage::Format_Grayscale8);
        if (!image.allGray()) {
            to_gray(image);
        }
        QImage grayImage = image.convertToFormat(QImage::Format_Grayscale8);
        toBinImage(grayImage);
        auto contours = detectContoursV(grayImage);
        std::vector<double> magnituteSpectrum = applyFFT(contours);

        uint8_t predictedClass = classifier.predict(magnituteSpectrum);
        guesses.set(i, testLabels[i] == predictedClass);
    }

    qDebug() << "correct" << guesses.count() << "out of" << testLabels.size();

    // shave off high frequencies and observe the change in accuracy

    qDebug() << "Cleaning up MNIST dataset files...";
    for (auto filePath : filePaths) {
        std::remove(filePath.c_str());
    }
    downloadDialog.close();
    // qDebug() << "elems" << QApplication::topLevelWidgets();
    app.quit();
    return 0;
}

#include "main.moc"
