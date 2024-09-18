#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_clean_clicked()
{
    ui->n->clear();
    ui->x0->clear();
    ui->xf->clear();
    ui->rhol->clear();
    ui->ul->clear();
    ui->pl->clear();
    ui->rhor->clear();
    ui->ur->clear();
    ui->pr->clear();
    ui->pd->clear();
    ui->gamma->clear();
    ui->tout->clear();
    ui->maxt->clear();
    ui->cfl->clear();
    ui->tol->clear();

    ui->problem->setCurrentIndex(-1);
    ui->BC1->setCurrentIndex(-1);
    ui->BC2->setCurrentIndex(-1);

    ui->radioButton->setChecked(true);
    ui->show->setChecked(true);

}

void MainWindow::on_problem_activated(int index)
{
    if(index==0)
    {
        ui->n->setText("100");
        ui->x0->setText("0.0");
        ui->xf->setText("1.0");
        ui->rhol->setText("1.0");
        ui->ul->setText("0.0");
        ui->pl->setText("1.0");
        ui->rhor->setText("0.125");
        ui->ur->setText("0.0");
        ui->pr->setText("0.1");
        ui->pd->setText("0.5");
        ui->gamma->setText("1.4");
        ui->tout->setText("0.25");
        ui->maxt->setText("100000");
        ui->cfl->setText("0.9");
        ui->tol->setText("0.000001");
        ui->BC1->setCurrentIndex(0);
        ui->BC2->setCurrentIndex(0);
        ui->show->setChecked(true);
    }
    if(index==1)
    {
        ui->n->setText("100");
        ui->x0->setText("0.0");
        ui->xf->setText("1.0");
        ui->rhol->setText("1.0");
        ui->ul->setText("-2.0");
        ui->pl->setText("0.4");
        ui->rhor->setText("1.0");
        ui->ur->setText("2.0");
        ui->pr->setText("0.4");
        ui->pd->setText("0.5");
        ui->gamma->setText("1.4");
        ui->tout->setText("0.15");
        ui->maxt->setText("100000");
        ui->cfl->setText("0.9");
        ui->tol->setText("0.000001");
        ui->BC1->setCurrentIndex(0);
        ui->BC2->setCurrentIndex(0);
        ui->show->setChecked(true);
    }
    if(index==2)
    {
        ui->n->setText("100");
        ui->x0->setText("0.0");
        ui->xf->setText("1.0");
        ui->rhol->setText("1.0");
        ui->ul->setText("-19.59745");
        ui->pl->setText("1000.0");
        ui->rhor->setText("1.0");
        ui->ur->setText("-19.59745");
        ui->pr->setText("0.01");
        ui->pd->setText("0.8");
        ui->gamma->setText("1.4");
        ui->tout->setText("0.012");
        ui->maxt->setText("100000");
        ui->cfl->setText("0.9");
        ui->tol->setText("0.000001");
        ui->BC1->setCurrentIndex(0);
        ui->BC2->setCurrentIndex(0);
        ui->show->setChecked(true);
    }
    if(index==3)
    {
        ui->n->setText("100");
        ui->x0->setText("0.0");
        ui->xf->setText("1.0");
        ui->rhol->setText("1.0");
        ui->ul->setText("0.0");
        ui->pl->setText("0.01");
        ui->rhor->setText("1.0");
        ui->ur->setText("0.0");
        ui->pr->setText("100.0");
        ui->pd->setText("0.5");
        ui->gamma->setText("1.4");
        ui->tout->setText("0.035");
        ui->maxt->setText("100000");
        ui->cfl->setText("0.9");
        ui->tol->setText("0.000001");
        ui->BC1->setCurrentIndex(0);
        ui->BC2->setCurrentIndex(0);
        ui->show->setChecked(true);
    }
    if(index==4)
    {
        ui->n->setText("100");
        ui->x0->setText("0.0");
        ui->xf->setText("1.0");
        ui->rhol->setText("5.99924");
        ui->ul->setText("19.5975");
        ui->pl->setText("460.894");
        ui->rhor->setText("5.99242");
        ui->ur->setText("-6.19633");
        ui->pr->setText("46.0950");
        ui->pd->setText("0.5");
        ui->gamma->setText("1.4");
        ui->tout->setText("0.035");
        ui->maxt->setText("100000");
        ui->cfl->setText("0.9");
        ui->tol->setText("0.000001");
        ui->BC1->setCurrentIndex(0);
        ui->BC2->setCurrentIndex(0);
        ui->show->setChecked(true);
    }
}

void MainWindow::on_calculate_clicked()
{
    int i[5], k, err;
    double t, et, d[14];

    i[0] = (ui->n->toPlainText()).toInt();
    i[1] = ui->BC1->currentIndex();
    i[2] = ui->BC2->currentIndex();
    if(ui->radioButton->isChecked())
        i[3] = 0;
    else if(ui->radioButton_2->isChecked())
        i[3] = 1;
    else if(ui->radioButton_3->isChecked())
        i[3] = 2;
    else if(ui->radioButton_4->isChecked())
        i[3] = 3;
    else if(ui->radioButton_5->isChecked())
        i[3] = 4;
    if(ui->show->isChecked()) i[4] = 1;
    else i[4] = 0;

    d[0] = (ui->rhol->toPlainText()).toDouble();
    d[1] = (ui->ul->toPlainText()).toDouble();
    d[2] = (ui->pl->toPlainText()).toDouble();
    d[3] = (ui->rhor->toPlainText()).toDouble();
    d[4] = (ui->ur->toPlainText()).toDouble();
    d[5] = (ui->pr->toPlainText()).toDouble();
    d[6] = (ui->x0->toPlainText()).toDouble();
    d[7] = (ui->xf->toPlainText()).toDouble();
    d[8] = (ui->pd->toPlainText()).toDouble();
    //if( (d[7] <= x1) || (d[7] >= x2) )
    d[9] = (ui->gamma->toPlainText()).toDouble();
    d[10] = (ui->tout->toPlainText()).toDouble();
    d[11] = (ui->maxt->toPlainText()).toDouble();
    d[12] = (ui->cfl->toPlainText()).toDouble();
    d[13] = (ui->tol->toPlainText()).toDouble();

    euler(d, i, &et, &t, &k, &err);

    if(i[3] == 1) ui->log->append("Roe Solver");
    else if(i[3] == 2) ui->log->append("HLL Solver");
    else if(i[3] == 3) ui->log->append("HLLC Solver");
    else if(i[3] == 4) ui->log->append("Rusanov Solver");
    else ui->log->append("Exact RP Solver");

    if(err == -1) ui->log->append("Error: Vacuum is generated \n "
                                  "           program was stopped\n");
    else
    {
        QString iter(QString("Iterations number: %1").arg(k));
        ui->log->append(iter);
        QString time(QString("Time:              %1 s").arg(t));
        ui->log->append(time);
        QString etime(QString("Elapsed time:     %1 s\n").arg(et));
        ui->log->append(etime);
    }

    if(err != -1) system("gnuplot gnuplot.gp && evince data.pdf");

}
