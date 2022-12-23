from PyQt5 import QtCore, QtWidgets
import pyqtgraph as pg
from math import log1p, floor
from statistics import variance, stdev, fmean, median, mode
from scipy.stats import norm, t, chi2



def reaction(object):
    object.calc_butt.clicked.connect(lambda: calc_sample(object))
    
    object.calc_gen_mean.clicked.connect(lambda: conf_calc_1(object))
    object.calc_gen_var.clicked.connect(lambda: conf_calc_2(object))
    
    object.calc_but.clicked.connect(lambda: gip_calc_1(object))
    object.calc_but_2.clicked.connect(lambda: gip_calc_2(object))
                                                     
def calc_sample(object):
    y_pos = 0
    
    #Создаем блок с ответами
    object.output_box = QtWidgets.QScrollArea(object.range_charact)
    object.output_box.setGeometry(QtCore.QRect(20, 335, 760, 200))
    object.output_box.setStyleSheet("background-color: rgb(234, 234, 234);")
    object.output_box.setWidgetResizable(True)
    object.output_box.setObjectName("output_box")
    object.outputBoxWidgetContents = QtWidgets.QWidget()
    object.outputBoxWidgetContents.setGeometry(QtCore.QRect(0, 0, 416, 196))
    object.outputBoxWidgetContents.setObjectName("outputBoxWidgetContents")
    object.gridLayout_2 = QtWidgets.QGridLayout(object.outputBoxWidgetContents)
    object.gridLayout_2.setObjectName("gridLayout_2")
    object.output_box.setWidget(object.outputBoxWidgetContents)
    
    #Мы отправляем текущую позицию по y. Потом задаем новое значение после выполнения функции как ее результат 
    
    #Дискретный вариационный ряд
    if object.discrete_variation_series.isChecked():
        y_pos = discrete_variation_series_f(object, y_pos)
    
    # Ряд распределения относительных частот
    if object.relative_frequency_distribution_series.isChecked():
        y_pos = relative_frequency_distribution_series_f(object, y_pos)
        
    # Полигон частот    
    if object.frequency_polygon.isChecked():
        y_pos = frequency_polygon_f(object, y_pos)
        
    # Эмпирическая функция распределния  
    if object.empirical_distribution_function.isChecked():
        y_pos = empirical_distribution_function_f(object, y_pos)
        
    # График эмпирической функции распределния  
    if object.plot_empirical_distribution_function.isChecked():
        y_pos = empirical_distribution_function_graph_f(object, y_pos)
    
    # Дисперсия  
    if object.variance.isChecked():
        y_pos = variance_f(object, y_pos)
        
    # Выборочное среднее 
    if object.sample_mean.isChecked():
        y_pos = sample_mean_f(object, y_pos)
        
    # Выборочное стандартное отклонение 
    if object.sample_standard_deviation.isChecked():
        y_pos = sample_standard_deviation_f(object, y_pos)
        
    # Исправленная дисперсия 
    if object.correct_variance.isChecked():
        y_pos = correct_variance_f(object, y_pos)
        
    # Мода 
    if object.moda.isChecked():
        y_pos = moda_f(object, y_pos)
        
    # Медиана  
    if object.mean.isChecked():
        y_pos = mean_f(object, y_pos)
        
    # Коэф вариаций 
    if object.var_k.isChecked():
        y_pos = var_k_f(object, y_pos)
        
    #Интервалы
        
    #Дискретный вариационный ряд
    if object.discrete_variation_series_interval.isChecked():
        y_pos = discrete_variation_series_f(object, y_pos)
    
    # Ряд распределения относительных частот
    if object.relative_frequency_distribution_series_interval.isChecked():
        y_pos = relative_frequency_distribution_series_interval_f(object, y_pos)
        
    # Полигон частот    
    if object.frequency_polygon_interval.isChecked():
        y_pos = frequency_polygon_interval_f(object, y_pos)
        
    # Эмпирическая функция распределния  
    if object.empirical_distribution_function_interval.isChecked():
        y_pos = empirical_distribution_function_f(object, y_pos)
        
    # График эмпирической функции распределния  
    if object.plot_empirical_distribution_function_interval.isChecked():
        y_pos = empirical_distribution_function_graph_f(object, y_pos)
    
    # Дисперсия  
    if object.dispersion_interval.isChecked():
        y_pos = variance_f(object, y_pos)
        
    # Выборочное среднее 
    if object.sample_mean_interval.isChecked():
        y_pos = sample_mean_f(object, y_pos)
        
    # Выборочное стандартное отклонение 
    if object.sample_standard_deviation_interval.isChecked():
        y_pos = sample_standard_deviation_f(object, y_pos)
        
    # Исправленная дисперсия 
    if object.correct_variance_interval.isChecked():
        y_pos = correct_variance_f(object, y_pos)
        
    # Мода 
    if object.moda_interval.isChecked():
        y_pos = moda_f(object, y_pos)
        
    # Медиана  
    if object.mean_interval.isChecked():
        y_pos = mean_f(object, y_pos)
        
    # Коэф вариаций 
    if object.var_k_interval.isChecked():
        y_pos = var_k_f(object, y_pos)
        

def conf_calc_1(object):
     
    conf = float(object.input_conf_gen_mean.toPlainText())
    n = float(object.input_n_gen_mean.toPlainText())
    v_mean = float(object.input_v_mean_gen_mean.toPlainText()) 
    
    if object.otkl_rbut.isChecked():
        otkl = float(object.input_corr_var_otkl_gen_mean.toPlainText())
        
        conf_z = min([-0.5 + norm.cdf(i/100) for i in range(1, 500)], key=lambda x: abs(x-conf/2))
        t_y = float([i/100 for i in range(1, 500) if -0.5 + norm.cdf(i/100) == conf_z][0])
        teta = (t_y*otkl)/(n**0.5)
        object.ensw_gen_mean.setText(f'Ответ: ({float("{:.3f}".format(v_mean - teta))};{(float("{:.3f}".format(v_mean+teta)))})')
        
    if object.corr_var_rbut.isChecked():
        corr_var = float(object.input_corr_var_otkl_gen_mean.toPlainText())
        t_y = t.ppf((1 +conf)/2, n-1)
        teta = (t_y*corr_var)/(n**0.5)
        object.ensw_gen_mean.setText(f'Ответ: ({float("{:.3f}".format(v_mean - teta))};{(float("{:.3f}".format(v_mean+teta)))})')
        
        
def conf_calc_2(object):
    conf = float(object.input_conf_gen_var.toPlainText())
    n = float(object.input_n_gen_var.toPlainText())
    s = float(object.input_corr_var_gen_var.toPlainText())
    
    x_1 = float("{:.3f}".format(chi2.ppf((1-conf)/2, n-1)))
    x_2 = float("{:.3f}".format(chi2.ppf((1+conf)/2, n-1)))
    
    print(x_1, x_2)
    
    f1 = ((n-1)*(s**2))/(x_2)
    f2 = ((n-1)*(s**2))/(x_1)
    f3 = (((n-1)**0.5)*s)/(x_2**0.5)
    f4 = (((n-1)**0.5)*s)/(x_1**0.5)
    print(f1, f2)
    
    if object.gen_var_rbut.isChecked():
        object.answ_gen_var.setText(f'Ответ: ({float("{:.3f}".format(f1))};{(float("{:.3f}".format(f2)))})')
        
    if object.stotkl_rbut_gen_var.isChecked():
        object.answ_gen_var.setText(f'Ответ: ({float("{:.3f}".format(f3))};{(float("{:.3f}".format(f4)))})')
        
        
def gip_calc_1(object):
    n = float(object.input_n_gip.toPlainText())
    ur_zn = float(object.input_ur_zn_gip.toPlainText())
    v_mean = float(object.input_sam_mean_gip.toPlainText())
    a = float(object.input_a_gip.toPlainText())
    gip0 = object.gip_lst.currentText()
    
    if object.gen_var_gip_rb.isChecked():
        gen_var = float(object.input_gen_var_corr_var_gip.toPlainText())
        if gip0=='H1: a < a₀':
            #Левосторонняя область
            #Критерий U, функция Лапласа
            f_zn = min([-0.5 + norm.cdf(i/100) for i in range(1, 500)], key=lambda x: abs(x-((1-2*ur_zn)/2)))
            f_cr = float([i/100 for i in range(1, 500) if -0.5 + norm.cdf(i/100) == f_zn][0])
            u = ((v_mean - a) * (n**0.5))/(gen_var**0.5)
            if -f_cr < u:
                object.gen_mean_answ.setText(f'Ответ: гипотеза на уровне значимости {ur_zn} принимается')
            elif -f_cr > u:
                object.gen_mean_answ.setText(f'Ответ: гипотеза на уровне значимости {ur_zn} отвергается')

        elif gip0=='H1: a > a₀':
            #правосторонняя область
            #Критерий U, функция Лапласа
            f_zn = min([-0.5 + norm.cdf(i/100) for i in range(1, 500)], key=lambda x: abs(x-((1-2*ur_zn)/2)))
            f_cr = float([i/100 for i in range(1, 500) if -0.5 + norm.cdf(i/100) == f_zn][0])
            u = ((v_mean - a) * (n**0.5))/(gen_var**0.5)
            if f_cr < u:
                object.gen_mean_answ.setText(f'Ответ: гипотеза на уровне значимости {ur_zn} отвергается')
            elif f_cr > u:
                object.gen_mean_answ.setText(f'Ответ: гипотеза на уровне значимости {ur_zn} принимается')
        else:
             #двусторонняя область
            #Критерий U, функция Лапласа
            f_zn = min([-0.5 + norm.cdf(i/100) for i in range(1, 500)], key=lambda x: abs(x-((1-ur_zn)/2)))
            f_cr = float([i/100 for i in range(1, 500) if -0.5 + norm.cdf(i/100) == f_zn][0])
            u = ((v_mean - a) * (n**0.5))/(gen_var**0.5)
            if -f_cr < u < f_cr:
                object.gen_mean_answ.setText(f'Ответ: гипотеза на уровне значимости {ur_zn} принимается')
            else:
                object.gen_mean_answ.setText(f'Ответ: гипотеза на уровне значимости {ur_zn} отвергается')
        
    elif object.corr_var_gip_rb.isChecked():
        corr_var = float(object.input_gen_var_corr_var_gip.toPlainText())
        
        # Левостороняя область
        if gip0=='H1: a < a₀':
            t_cr = t.ppf(1-ur_zn, n-1)
            t_n = ((v_mean - a) * (n**0.5))/(corr_var**0.5)
            if -t_cr < t_n:
                object.gen_mean_answ.setText(f'Ответ: гипотеза на уровне значимости {ur_zn} принимается')
            elif -t_cr > t_n:
                object.gen_mean_answ.setText(f'Ответ: гипотеза на уровне значимости {ur_zn} отвергается')  
            
        # правостороняя область   
        elif gip0=='H1: a > a₀':
            t_cr = t.ppf(1-ur_zn, n-1)
            t_n = ((v_mean - a) * (n**0.5))/(corr_var**0.5)
            if t_cr < t_n:
                object.gen_mean_answ.setText(f'Ответ: гипотеза на уровне значимости {ur_zn} отвергается')
            elif t_cr > t_n:
                object.gen_mean_answ.setText(f'Ответ: гипотеза на уровне значимости {ur_zn} принимается')
        # двусторонняя область
        else:
            t_cr = t.ppf(1-ur_zn/2, n-1)
            t_n = ((v_mean - a) * (n**0.5))/(corr_var**0.5)
            if -t_cr < t_n < t_cr:
                object.gen_mean_answ.setText(f'Ответ: гипотеза на уровне значимости {ur_zn} принимается')
            else:
                object.gen_mean_answ.setText(f'Ответ: гипотеза на уровне значимости {ur_zn} отвергается')
    
    
def gip_calc_2(object):
    n = float(object. input_n_gip_2.toPlainText())
    sam_mean_1 = float(object.input_n_gip_5.toPlainText())
    m = float(object.input_n_gip_3.toPlainText())
    sam_mean_2 = float(object.input_n_gip_4.toPlainText())
    var_1 = float(object.input_n_gip_6.toPlainText())
    var_2 = float(object.input_n_gip_7.toPlainText())
    ur_zn = float(object.input_n_gip_8.toPlainText())
    
    gip0 = object.gip_lst_2.currentText()
    
    if gip0=='H1: Ẋв₁ < Ẋв₂':
        #Левосторонняя область
        #Критерий U, функция Лапласа
        f_zn = min([-0.5 + norm.cdf(i/100) for i in range(1, 500)], key=lambda x: abs(x-((1-2*ur_zn)/2)))
        f_cr = float([i/100 for i in range(1, 500) if -0.5 + norm.cdf(i/100) == f_zn][0])
        z_n = (sam_mean_1 - sam_mean_2)/(((var_1/n) + (var_2/m) )**0.5)
        if -f_cr < z_n:
            object.gen_mean_answ_2.setText(f'Ответ: гипотеза на уровне значимости {ur_zn} принимается')
        elif -f_cr > z_n:
            object.gen_mean_answ_2.setText(f'Ответ: гипотеза на уровне значимости {ur_zn} отвергается')

    elif gip0=='H1: Ẋв₁ > Ẋв₂ ':
        #правосторонняя область
        #Критерий U, функция Лапласа
        f_zn = min([-0.5 + norm.cdf(i/100) for i in range(1, 500)], key=lambda x: abs(x-((1-2*ur_zn)/2)))
        f_cr = float([i/100 for i in range(1, 500) if -0.5 + norm.cdf(i/100) == f_zn][0])
        z_n = (sam_mean_1 - sam_mean_2)/(((var_1/n) + (var_2/m) )**0.5)
        if f_cr < z_n:
            object.gen_mean_answ_2.setText(f'Ответ: гипотеза на уровне значимости {ur_zn} отвергается')
        elif f_cr > z_n:
            object.gen_mean_answ_2.setText(f'Ответ: гипотеза на уровне значимости {ur_zn} принимается')
    else:
         #двусторонняя область
        #Критерий U, функция Лапласа
        f_zn = min([-0.5 + norm.cdf(i/100) for i in range(1, 500)], key=lambda x: abs(x-((1-ur_zn)/2)))
        f_cr = float([i/100 for i in range(1, 500) if -0.5 + norm.cdf(i/100) == f_zn][0])
        z_n = (sam_mean_1 - sam_mean_2)/(((var_1/n) + (var_2/m) )**0.5)
        if -f_cr < z_n < f_cr:
            object.gen_mean_answ_2.setText(f'Ответ: гипотеза на уровне значимости {ur_zn} принимается')
        else:
            object.gen_mean_answ_2.setText(f'Ответ: гипотеза на уровне значимости {ur_zn} отвергается')
    
    

    
    
    
def discrete_variation_series_f(object, y_pos):
    sample = object.input_box.toPlainText().split(' ')
    dict_sample = {}
    col1 = list(set(sample))
    for k in col1:
        dict_sample[k] = sample.count(k)
    
    object.title_dis = QtWidgets.QLabel(object.outputBoxWidgetContents)
    object.title_dis.setStyleSheet("color: rgb(0, 0, 0);\n"
"background-color: rgb(234, 234, 234);")
    object.title_dis.setObjectName("title_dis")
    object.title_dis.setGeometry(QtCore.QRect(20, y_pos+20, 800, 15))
    object.title_dis.setText("Дискретный вариационный ряд:")
    object.gridLayout_2.addWidget(object.title_dis)
    y_pos += 35
    
    object.dis_var_r_table = QtWidgets.QTableWidget(object.outputBoxWidgetContents)
    object.dis_var_r_table.setObjectName("dis_var_r_table")
    object.dis_var_r_table.setGeometry(QtCore.QRect(20, y_pos+20, 800, 60))
    object.dis_var_r_table.setColumnCount(len(col1))
    object.dis_var_r_table.setRowCount(2)
    object.dis_var_r_table.setVerticalHeaderLabels(['xᵢ', 'nᵢ'])
    r = 0
    for par in sorted(dict_sample.items(), key=lambda x: x[0]):
        object.dis_var_r_table.setItem(0, r, QtWidgets.QTableWidgetItem(par[0]))
        object.dis_var_r_table.setItem(1, r, QtWidgets.QTableWidgetItem(par[1]))
        r+=1
    object.gridLayout_2.addWidget(object.dis_var_r_table)
    object.output_box.show()
    
    y_pos += 80
    return y_pos


def relative_frequency_distribution_series_f(object, y_pos):
    sample = object.input_box.toPlainText().split(' ')
    dict_sample = {}
    col1 = list(set(sample))
    for k in col1:
        dict_sample[k] = sample.count(k)
    
    object.title_frequency = QtWidgets.QLabel(object.outputBoxWidgetContents)
    object.title_frequency.setStyleSheet("color: rgb(0, 0, 0);\n"
"background-color: rgb(234, 234, 234);")
    object.title_frequency.setObjectName("title_frequency")
    object.title_frequency.setGeometry(QtCore.QRect(20, y_pos+20, 800, 15))
    object.title_frequency.setText("Ряд распределения относительных частот:")
    object.gridLayout_2.addWidget(object.title_frequency)
    
    y_pos += 35
    
    object.frequency_table = QtWidgets.QTableWidget(object.outputBoxWidgetContents)
    object.frequency_table.setObjectName("frequency_table")
    object.frequency_table.setGeometry(QtCore.QRect(20, y_pos+20, 800, 60))
    object.frequency_table.setColumnCount(len(col1))
    object.frequency_table.setRowCount(2)
    object.frequency_table.setVerticalHeaderLabels(['xᵢ', 'wᵢ'])
    r = 0
    for par in sorted(dict_sample.items(), key=lambda x: x[0]):
        object.frequency_table.setItem(0, r, QtWidgets.QTableWidgetItem(par[0]))
        object.frequency_table.setItem(1, r, QtWidgets.QTableWidgetItem(str(float('{:.3f}'.format(par[1] / len(sample))))))
        r+=1
    object.gridLayout_2.addWidget(object.frequency_table)
    object.output_box.show()
    
    y_pos += 80
    
    return y_pos
    
    
def frequency_polygon_f(object, y_pos):
    sample = object.input_box.toPlainText().split(' ')
    dict_sample = {}
    col1 = list(set(sample))
    for k in col1:
        dict_sample[k] = sample.count(k)
    
    object.frequency_polygon_graph = pg.PlotWidget(object.outputBoxWidgetContents)
    object.frequency_polygon_graph.setStyleSheet("color: rgb(0, 0, 0);\n"
"background-color: rgb(234, 234, 234);")
    object.frequency_polygon_graph.setObjectName("frequency_polygon_graph")
    object.frequency_polygon_graph.setGeometry(QtCore.QRect(20, y_pos+20, 600, 120))
    object.frequency_polygon_graph.setFixedSize(600, 120)
    values = [float(par[0]) for par in sorted(dict_sample.items(), key=lambda x: x[0])]
    n = [float(par[1]/len(sample)) for par in sorted(dict_sample.items(), key=lambda x: x[0])]
    
    object.frequency_polygon_graph.setBackground((234, 234, 234))
    object.frequency_polygon_graph.setTitle("Полигон частот", size="10pt", color='black')
    styles = {"color": "#f00", "font-size": "10px"}
    object.frequency_polygon_graph.setLabel("left", "wᵢ", **styles)
    object.frequency_polygon_graph.setLabel("bottom", "'xᵢ", **styles)
    object.frequency_polygon_graph.addLegend()
    object.frequency_polygon_graph.showGrid(x=True, y=True)
    pen = pg.mkPen(color=(46, 139, 87))
    
    object.frequency_polygon_graph.plot(values, n,  pen=pen)
    object.gridLayout_2.addWidget(object.frequency_polygon_graph)
    
    y_pos += 120
    object.output_box.show()
    return y_pos
    
    
def empirical_distribution_function_f(object, y_pos):
    object.title_emp = QtWidgets.QLabel(object.outputBoxWidgetContents)
    object.title_emp.setStyleSheet("color: rgb(0, 0, 0);\n"
"background-color: rgb(234, 234, 234);")
    object.title_emp.setObjectName("title_emp")
    object.title_emp.setGeometry(QtCore.QRect(20, y_pos+20, 800, 15))
    object.title_emp.setText("Эмпирическая функция распределения")
    object.gridLayout_2.addWidget(object.title_emp)
    y_pos += 35
    
    sample = sorted(object.input_box.toPlainText().split(' '))
    dict_sample = {}
    col1 = list(set(sample))
    for k in col1:
        dict_sample[k] = sample.count(k)
        
    k = floor(1 + 3.322*log1p(len(sample)))
    sam_min, sam_max = float(min(sample)), float(max(sample))
    r = sam_max - sam_min
    h = r/k
    lst_interv = [[sam_min + i*h, sam_min + (i+1)*h] for i in range(k)]
    for interval in lst_interv:
        interval.append(len([el for el in sample if float(el) >= interval[0] and float(el) < interval[1]])/len(sample))
    lst_interv[-1][2] += len([el for el in sample if float(el) == lst_interv[-1][1]])/len(sample)
    
    for interval in range(len(lst_interv)):
        if interval > 0:
            lst_interv[interval][2] += lst_interv[interval-1][2]
    
    object.emp_table = QtWidgets.QTableWidget(object.outputBoxWidgetContents)
    object.emp_table.setObjectName("emp_table")
    object.emp_table.setGeometry(QtCore.QRect(20, y_pos+20, 800, 100))
    object.emp_table.setFixedSize(800, 100)
    object.emp_table.setColumnCount(2)
    object.emp_table.setRowCount(len(lst_interv))
    object.emp_table.setHorizontalHeaderLabels(['При', 'F(x)'])
    
    r = 0
    for interval in lst_interv:
        object.emp_table.setItem(r, 0, QtWidgets.QTableWidgetItem(f'Х от {float("{:.3f}".format(interval[0]))} до {float("{:.3f}".format(interval[1]))}'))
        object.emp_table.setItem(r, 1, QtWidgets.QTableWidgetItem(f'{interval[2]}'))
        r+=1
    object.gridLayout_2.addWidget(object.emp_table)
    object.output_box.show()
    
    y_pos += 120
    object.output_box.show()
    
    return y_pos

    
def empirical_distribution_function_graph_f(object, y_pos):
    sample = sorted(object.input_box.toPlainText().split(' '))
   
    k = floor(1 + 3.322*log1p(len(sample)))
    sam_min, sam_max = float(min(sample)), float(max(sample))
    r = sam_max - sam_min
    h = r/k
    lst_interv = [[sam_min + i*h, sam_min + (i+1)*h] for i in range(k)]
    for interval in lst_interv:
        interval.append(len([el for el in sample if float(el) >= interval[0] and float(el) < interval[1]])/len(sample))
    lst_interv[-1][2] += len([el for el in sample if float(el) == lst_interv[-1][1]])/len(sample)
    
    for interval in range(len(lst_interv)):
        if interval > 0:
            lst_interv[interval][2] += lst_interv[interval-1][2]
            
    
    object.emp_graph = pg.PlotWidget(object.outputBoxWidgetContents)
    object.emp_graph.setStyleSheet("color: rgb(0, 0, 0);\n"
"background-color: rgb(234, 234, 234);")
    object.emp_graph.setObjectName("emp_graph")
    object.emp_graph.setGeometry(QtCore.QRect(20, y_pos+20, 400, 100))
    object.emp_graph.setFixedSize(600, 120)
    object.emp_graph.setBackground((234, 234, 234))
    object.emp_graph.setTitle("График эмпирической функции распределения", size="10pt", color='black')
    styles = {"color": "#f00", "font-size": "10px"}
    object.emp_graph.setLabel("left", "Относ. частота (wᵢ)", **styles)
    object.emp_graph.setLabel("bottom", "Диапазон (x)", **styles)
    object.emp_graph.addLegend()
    object.emp_graph.showGrid(x=True, y=True)
    pen = pg.mkPen(color=(46, 139, 87))
        
    for interval in lst_interv:
        object.emp_graph.plot([interval[0], interval[1]], [interval[2], interval[2]],  pen=pen)
    object.gridLayout_2.addWidget(object.emp_graph)
    
    y_pos += 120
    object.output_box.show()
    
    return y_pos

def variance_f(object, y_pos):
    object.var_answ = QtWidgets.QLabel(object.outputBoxWidgetContents)
    object.var_answ.setStyleSheet("color: rgb(0, 0, 0);\n"
"background-color: rgb(234, 234, 234);")
    object.var_answ.setObjectName("var_answ")
    object.var_answ.setGeometry(QtCore.QRect(20, y_pos+20, 800, 15))
    object.var_answ.setText(f"Дисперсия: {variance([float(el) for el in sorted(object.input_box.toPlainText().split(' '))])}")
    object.gridLayout_2.addWidget(object.var_answ)
    y_pos += 35
    
    object.output_box.show()
    
    return y_pos


def sample_mean_f(object, y_pos):
    object.sample_mean_answ = QtWidgets.QLabel(object.outputBoxWidgetContents)
    object.sample_mean_answ.setStyleSheet("color: rgb(0, 0, 0);\n"
"background-color: rgb(234, 234, 234);")
    object.sample_mean_answ.setObjectName("sample_mean_answ")
    object.sample_mean_answ.setGeometry(QtCore.QRect(20, y_pos+20, 800, 15))
    object.sample_mean_answ.setText(f"Выборочное среднее: {fmean([float(el) for el in sorted(object.input_box.toPlainText().split(' '))])}")
    object.gridLayout_2.addWidget(object.sample_mean_answ)
    y_pos += 35
    
    object.output_box.show()
    
    return y_pos

def sample_standard_deviation_f(object, y_pos):
    object.sample_standard_deviation_answ = QtWidgets.QLabel(object.outputBoxWidgetContents)
    object.sample_standard_deviation_answ.setStyleSheet("color: rgb(0, 0, 0);\n"
"background-color: rgb(234, 234, 234);")
    object.sample_standard_deviation_answ.setObjectName("sample_standard_deviation_answ")
    object.sample_standard_deviation_answ.setGeometry(QtCore.QRect(20, y_pos+20, 800, 15))
    object.sample_standard_deviation_answ.setText(f"Выборочное стандартное отклонение: {stdev([float(el) for el in sorted(object.input_box.toPlainText().split(' '))])}")
    object.gridLayout_2.addWidget(object.sample_standard_deviation_answ)
    y_pos += 35
    
    object.output_box.show()
    
    return y_pos


def correct_variance_f(object, y_pos):
    sample = sorted(object.input_box.toPlainText().split(' '))
    object.correct_variance_answ = QtWidgets.QLabel(object.outputBoxWidgetContents)
    object.correct_variance_answ.setStyleSheet("color: rgb(0, 0, 0);\n"
"background-color: rgb(234, 234, 234);")
    object.correct_variance_answ.setObjectName("correct_variance_answ")
    object.correct_variance_answ.setGeometry(QtCore.QRect(20, y_pos+20, 800, 15))
    object.correct_variance_answ.setText(f"Исправленная дисперсия: {(variance([float(el) for el in sample]) * (len(sample)/(len(sample) - 1)))**0.5}")
    object.gridLayout_2.addWidget(object.correct_variance_answ)
    y_pos += 35
    
    object.output_box.show()
    
    return y_pos


def moda_f(object, y_pos):
    object.moda_answ = QtWidgets.QLabel(object.outputBoxWidgetContents)
    object.moda_answ.setStyleSheet("color: rgb(0, 0, 0);\n"
"background-color: rgb(234, 234, 234);")
    object.moda_answ.setObjectName("moda_answ")
    object.moda_answ.setGeometry(QtCore.QRect(20, y_pos+20, 800, 15))
    object.moda_answ.setText(f"Мода: {mode([float(el) for el in sorted(object.input_box.toPlainText().split(' '))])}")
    object.gridLayout_2.addWidget(object.moda_answ)
    y_pos += 35
    
    object.output_box.show()
    
    return y_pos


def mean_f(object, y_pos):
    object.mean_answ = QtWidgets.QLabel(object.outputBoxWidgetContents)
    object.mean_answ.setStyleSheet("color: rgb(0, 0, 0);\n"
"background-color: rgb(234, 234, 234);")
    object.mean_answ.setObjectName("mean_answ")
    object.mean_answ.setGeometry(QtCore.QRect(20, y_pos+20, 800, 15))
    object.mean_answ.setText(f"Медиана: {median([float(el) for el in sorted(object.input_box.toPlainText().split(' '))])}")
    object.gridLayout_2.addWidget(object.mean_answ)
    y_pos += 35
    
    object.output_box.show()
    
    return y_pos


def var_k_f(object, y_pos):
    object.var_k_answ = QtWidgets.QLabel(object.outputBoxWidgetContents)
    object.var_k_answ.setStyleSheet("color: rgb(0, 0, 0);\n"
"background-color: rgb(234, 234, 234);")
    object.var_k_answ.setObjectName("var_k_answ")
    object.var_k_answ.setGeometry(QtCore.QRect(20, y_pos+20, 800, 15))
    object.var_k_answ.setText(f"Коэффициент вариаций: {(stdev([float(el) for el in sorted(object.input_box.toPlainText().split(' '))])/fmean([float(el) for el in sorted(object.input_box.toPlainText().split(' '))]))*100} %")
    object.gridLayout_2.addWidget(object.var_k_answ)
    y_pos += 35
    
    object.output_box.show()
    
    return y_pos

    
def  relative_frequency_distribution_series_interval_f(object, y_pos):
    sample = object.input_box.toPlainText().split(' ')

    object.title_frequency_int = QtWidgets.QLabel(object.outputBoxWidgetContents)
    object.title_frequency_int.setStyleSheet("color: rgb(0, 0, 0);\n"
"background-color: rgb(234, 234, 234);")
    object.title_frequency_int.setObjectName("title_frequency_int")
    object.title_frequency_int.setGeometry(QtCore.QRect(20, y_pos+20, 800, 15))
    object.title_frequency_int.setText("Интервальный ряд распределения относительных частот:")
    object.gridLayout_2.addWidget(object.title_frequency_int)
    
    y_pos += 35
    
     # Вычислим, сколько должно быть интервалов
    k = floor(1 + 3.322*log1p(len(sample)))
    sam_min, sam_max = float(min(sample)), float(max(sample))
    r = sam_max - sam_min
    h = r/k
    lst_interv = [[sam_min + i*h, sam_min + (i+1)*h] for i in range(k)]
    for interval in lst_interv:
        interval.append(len([el for el in sample if float(el) >= interval[0] and float(el) < interval[1]])/len(sample))
    lst_interv[-1][2] += len([el for el in sample if float(el) == lst_interv[-1][1]])/len(sample)
    
    for interval in range(len(lst_interv)):
        if interval > 0:
            lst_interv[interval][2] += lst_interv[interval-1][2]
    
    object.frequency_table_int = QtWidgets.QTableWidget(object.outputBoxWidgetContents)
    object.frequency_table_int.setObjectName("frequency_table_int")
    object.frequency_table_int.setGeometry(QtCore.QRect(20, y_pos+20, 800, 60))
    object.frequency_table_int.setColumnCount(k)
    object.frequency_table_int.setRowCount(2)
    object.frequency_table_int.setVerticalHeaderLabels(['Интервалы', 'wᵢ'])
    r = 0
    for interval in lst_interv:
        object.frequency_table_int.setItem(0, r, QtWidgets.QTableWidgetItem(f'{interval[0]} - {interval[1]}'))
        object.frequency_table_int.setItem(1, r, QtWidgets.QTableWidgetItem(str(float('{:.3f}'.format(interval[2] / len(sample))))))
        r+=1
    object.gridLayout_2.addWidget(object.frequency_table_int)
    object.output_box.show()
    
    y_pos += 80
    
    return y_pos


def frequency_polygon_interval_f(object, y_pos):
    object.frequency_polygon_graph_int = pg.PlotWidget(object.outputBoxWidgetContents)
    object.frequency_polygon_graph_int.setStyleSheet("color: rgb(0, 0, 0);\n"
"background-color: rgb(234, 234, 234);")
    object.frequency_polygon_graph_int.setObjectName("frequency_polygon_graph_int")
    object.frequency_polygon_graph_int.setGeometry(QtCore.QRect(20, y_pos+20, 600, 120))
    object.frequency_polygon_graph_int.setFixedSize(600, 120)
    object.frequency_polygon_graph_int.setBackground((234, 234, 234))
    object.frequency_polygon_graph_int.setTitle("Полигон частот интервального ряда", size="10pt", color='black')
    styles = {"color": "#f00", "font-size": "10px"}
    object.frequency_polygon_graph_int.setLabel("left", "wᵢ", **styles)
    object.frequency_polygon_graph_int.setLabel("bottom", "'xᵢ", **styles)
    object.frequency_polygon_graph_int.addLegend()
    object.frequency_polygon_graph_int.showGrid(x=True, y=True)
    pen = pg.mkPen(color=(46, 139, 87))
    
    sample = sorted(object.input_box.toPlainText().split(' '))
    # Вычислим, сколько должно быть интервалов
    k = floor(1 + 3.322*log1p(len(sample)))
    sam_min, sam_max = float(min(sample)), float(max(sample))
    r = sam_max - sam_min
    h = r/k
    lst_interv = [[sam_min + i*h, sam_min + (i+1)*h] for i in range(k)]
    for interval in lst_interv:
        interval.append(len([el for el in sample if float(el) >= interval[0] and float(el) < interval[1]])/len(sample))
    lst_interv[-1][2] += len([el for el in sample if float(el) == lst_interv[-1][1]])/len(sample)
    
    object.frequency_polygon_graph_int.plot([sam_min + i*h for i in range(k+1)], [0] + [interval[2] for interval in lst_interv],  pen=pen)
    object.gridLayout_2.addWidget(object.frequency_polygon_graph_int)
    
    y_pos += 120
    object.output_box.show()
    return y_pos