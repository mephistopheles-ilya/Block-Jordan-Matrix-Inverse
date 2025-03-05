#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QAction>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QMessageBox>

#include "window.hpp"

#include <fenv.h>

int main (int argc, char *argv[]) {

    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);

    QApplication app (argc, argv);

    QMainWindow *window = new QMainWindow;
    QMenuBar *tool_bar = new QMenuBar (window);
    Window *graph_area = new Window (window);
    QAction *action;

    if (graph_area->parse_command_line (argc, argv)) {
        QMessageBox::warning (0, "Wrong input arguments!", 
                "Wrong input arguments!");
        return -1;
    }

    action = tool_bar->addAction ("&Change function", graph_area, SLOT (change_func ()));
    action->setShortcut (QString ("Ctrl+C"));

    action = tool_bar->addAction ("E&xit", window, SLOT (close ()));
    action->setShortcut (QString ("Ctrl+X"));

    tool_bar->setMaximumHeight (30);

    window->setMenuBar (tool_bar);
    window->setCentralWidget (graph_area);
    window->setWindowTitle ("Graph");

    window->show ();
    app.exec ();
    delete window;
    return 0;
}

#if 0
#include <QMainWindow>
#include <QPainter>
#include <QMenuBar>
#include <QApplication>
#include <QStyleFactory>
#include <cmath>

class MainWindow: public QMainWindow {
    Q_OBJECT
private:
        QString currentFunction;
        int numPoints;
        double xMin, xMax;
public:
        explicit MainWindow(QWidget* parent = nullptr): QMainWindow(parent), currentFunction("sin"), numPoints(1000), xMin(-10.), xMax(10.) {
            setWindowTitle("Math Function Plotter");
            resize(800, 600); // Set initial window size

            // Create a menu bar
            QMenuBar *menuBar = new QMenuBar(this);
            setMenuBar(menuBar);

            // Create a "Functions" menu
            QMenu *functionsMenu = menuBar->addMenu("Functions");

            // Add actions to the menu with icons and shortcuts
            QAction *sinAction = functionsMenu->addAction(QIcon::fromTheme("sin"), "Sine");
            sinAction->setShortcut(QKeySequence("Ctrl+S")); // Shortcut for Sine
            QAction *cosAction = functionsMenu->addAction(QIcon::fromTheme("cos"), "Cosine");
            cosAction->setShortcut(QKeySequence("Ctrl+C")); // Shortcut for Cosine
            QAction *quadraticAction = functionsMenu->addAction(QIcon::fromTheme("quadratic"), "Quadratic");
            quadraticAction->setShortcut(QKeySequence("Ctrl+Q")); // Shortcut for Quadratic

            // Add a separator
            functionsMenu->addSeparator();

            // Add an exit action
            QAction *exitAction = functionsMenu->addAction(QIcon::fromTheme("exit"), "Exit");
            exitAction->setShortcut(QKeySequence("Ctrl+E")); // Shortcut for Exit
            connect(exitAction, &QAction::triggered, this, &QMainWindow::close);

            // Create a "Settings" menu
            QMenu *settingsMenu = menuBar->addMenu("Settings");

            // Add actions to adjust the number of points
            QAction *increasePointsAction = settingsMenu->addAction("Increase Points");
            increasePointsAction->setShortcut(QKeySequence("Ctrl++")); // Shortcut for Increase Points
            QAction *decreasePointsAction = settingsMenu->addAction("Decrease Points");
            decreasePointsAction->setShortcut(QKeySequence("Ctrl+-")); // Shortcut for Decrease Points

            // Add actions to adjust the axes length
            QAction *increaseAxesAction = settingsMenu->addAction("Increase Axes Length");
            increaseAxesAction->setShortcut(QKeySequence("Ctrl+Shift++")); // Shortcut for Increase Axes Length
            QAction *decreaseAxesAction = settingsMenu->addAction("Decrease Axes Length");
            decreaseAxesAction->setShortcut(QKeySequence("Ctrl+Shift+-")); // Shortcut for Decrease Axes Length

            // Connect actions to slots
            connect(sinAction, &QAction::triggered, [this]() {
                    setFunction("sin");
                    });
            connect(cosAction, &QAction::triggered, [this]() {
                    setFunction("cos");
                    });
            connect(quadraticAction, &QAction::triggered, [this]() {
                    setFunction("quadratic");
                    });
            connect(increasePointsAction, &QAction::triggered, [this]() {
                    numPoints += 500; // Increase number of points
                    update();
                    });
            connect(decreasePointsAction, &QAction::triggered, [this]() {
                    numPoints = std::max(100, numPoints - 500); // Decrease number of points (minimum 100)
                    update();
                    });
            connect(increaseAxesAction, &QAction::triggered, [this]() {
                    xMin -= 5.0; // Increase X-axis range
                    xMax += 5.0;
                    update();
                    });
            connect(decreaseAxesAction, &QAction::triggered, [this]() {
                    xMin += 5.0; // Decrease X-axis range
                    xMax -= 5.0;
                    xMin = std::max(-100.0, xMin); // Prevent negative range
                    xMax = std::min(100.0, xMax);
                    update();
                    });
        }

        void setFunction(const QString &function) {
            currentFunction = function;
            update();
        }

        void paintEvent(QPaintEvent* ) override {
            QPainter painter(this);
            // Define the dimensions of the plot area
            double width = this->width();
            double height = this->height();
            double margin = 50; // Margin around the plot

            // Draw the axes
            painter.setPen(Qt::black);

            // X-axis
            painter.drawLine(margin, height / 2, width - margin, height / 2); // X-axis line
                                                                              // Draw arrow at the end of the X-axis
            painter.drawLine(width - margin, height / 2, width - margin - 10, height / 2 - 5); // Arrow part 1
            painter.drawLine(width - margin, height / 2, width - margin - 10, height / 2 + 5); // Arrow part 2

            // Y-axis
            painter.drawLine(width / 2, margin, width / 2, height - margin); // Y-axis line
                                                                             // Draw arrow at the end of the Y-axis
            painter.drawLine(width / 2, margin, width / 2 - 5, margin + 10); // Arrow part 1
            painter.drawLine(width / 2, margin, width / 2 + 5, margin + 10); // Arrow part 2

            // Define the Y range based on the selected function
            double yMin = -10.0, yMax = 10.0; // Default for sine/cosine
            if (currentFunction == "quadratic") {
                yMin = 0.0;
                yMax = xMax * xMax; // y = x^2 ranges from 0 to xMax^2
            }

            // Draw ticks and numbers on the X-axis
            double xStep = (xMax - xMin) / 10.0; // Step between ticks (10 ticks total)
            for (double x = xMin; x <= xMax; x += xStep) {
                double widgetX = margin + ((x - xMin) / (xMax - xMin)) * (width - 2 * margin);
                painter.drawLine(widgetX, height / 2 - 5, widgetX, height / 2 + 5); // Tick
                painter.drawText(widgetX - 10, height / 2 + 20, QString::number(x, 'f', 1)); // Number with 1 decimal place
            }

            // Draw ticks and numbers on the Y-axis
            double yStep = (yMax - yMin) / 10.0; // Step between ticks (10 ticks total)
            for (double y = yMin; y <= yMax; y += yStep) {
                double widgetY = height / 2 - ((y - yMin) / (yMax - yMin)) * (height - 2 * margin);
                painter.drawLine(width / 2 - 5, widgetY, width / 2 + 5, widgetY); // Tick
                painter.drawText(width / 2 + 10, widgetY + 5, QString::number(y, 'f', 1)); // Number with 1 decimal place
            }

            // Draw the selected function
            painter.setPen(Qt::blue);
            QPointF prevPoint;
            bool firstPoint = true;

            for (int i = 0; i < numPoints; ++i) {
                double x = xMin + (xMax - xMin) * i / (numPoints - 1);
                double y = 0.0;

                // Calculate y based on the selected function
                if (currentFunction == "sin") {
                    y = sin(x);
                } else if (currentFunction == "cos") {
                    y = cos(x);
                } else if (currentFunction == "quadratic") {
                    y = x * x; // y = x^2
                }

                // Map function coordinates to widget coordinates
                double widgetX = margin + ((x - xMin) / (xMax - xMin)) * (width - 2 * margin);
                double widgetY = ((y - yMin) / (yMax - yMin)) * (height - 2 * margin);

                QPointF currentPoint(widgetX, widgetY);

                if (!firstPoint) {
                    painter.drawLine(prevPoint, currentPoint); // Draw line between points
                }

                prevPoint = currentPoint;
                firstPoint = false;
            }
        }
        ~MainWindow() {};
};

#include <QApplication>

#include "window.hpp"

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);

    MainWindow window;
    window.show();

    return app.exec();
}
#endif
