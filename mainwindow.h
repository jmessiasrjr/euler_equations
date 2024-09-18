#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QObject>
#include <QLocale>
#include "euler.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void on_clean_clicked();

    void on_problem_activated(int index);

    void on_calculate_clicked();

private:
    Ui::MainWindow *ui;

};

#endif // MAINWINDOW_H
