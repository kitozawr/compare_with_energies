import sys, os
import numpy as np
from scipy.optimize import curve_fit
import scipy.stats as sps
import matplotlib.pyplot as plt
import shutil as sh
import argparse

#sys.path.append("./data_proc/data_proc_basics")
import data_proc as dp

def fit_function(x, a, b):
	return a * x + b

"""
def r_2(x, y, a, b):
	return 1 - np.sum((y - fit_function(x, a, b)) ** 2) / np.sum((y - y.mean()) ** 2)


def set_bounds(foldername_ac_splitted, inv, bounds_inv, bounds_not_inv):

	'''
	Select bounds for inverted or non-inverted acoustic signal.
	'''

	if not inv or 'kalibr' in foldername_ac_splitted or 'calibr' in foldername_ac_splitted:
		return bounds_not_inv
	return bounds_inv
"""

"""
def fit(fit_function, energies_parrots, maxima_new, bounds = (np.array([1e-5, -np.inf]), np.array([np.inf, np.inf]))):

	'''Get fitting parameters and r^2.'''

	popt, pcov = curve_fit(fit_function, energies_parrots, maxima_new, bounds=bounds)

	a = popt[0]
	b = popt[1]
	cur_r_2 = r_2(energies_parrots, maxima_new, a, b)
	print("Fit function: a = {}, b = {}, r_2 = {}".format(a,b,cur_r_2))
	a1, b1, rvalue = sps.linregress(energies_parrots, maxima_new)[0:3]
	print("Fit function linregress: a1 = {}, b1 = {}, rvalue^2 = {}".format(a1,b1,rvalue**2))
	print("Количество сопоставленных данных: {}".format(len(energies_parrots)))

	return(cur_r_2, a, b, a1, b1)


def fit(fit_function, energies_parrots, maxima_new, bounds = (np.array([1e-5, -np.inf]), np.array([np.inf, np.inf]))):

	'''Get fitting parameters and r^2.'''

	a, b, cur_r_2 = sps.linregress(energies_parrots, maxima_new)[0:3]
	print("Fit function: a = {}, b = {}, r_2 = {}".format(a,b,cur_r_2))
	a1, b1, rvalue = sps.linregress(energies_parrots, maxima_new)[0:3]
	print("Fit function linregress: a1 = {}, b1 = {}, rvalue^2 = {}".format(a1,b1,rvalue**2))
	print("Количество сопоставленных данных: {}".format(len(energies_parrots)))

	return(cur_r_2, a, b, a1, b1)

#Старая версия с polyfit.

def fit_and_check_shift(fit_function, energies_parrots, maxima_new, best_shift, shift, max_r_2, num_matched, set_limit, bounds = (np.array([1e-5, -np.inf]), np.array([np.inf, np.inf]))):

	popt, pcov = curve_fit(fit_function, energies_parrots, maxima_new, bounds=bounds)

	a = popt[0]
	b = popt[1]
	cur_r_2 = r_2(energies_parrots, maxima_new, a, b)
	print("Fit and check function: a = {}, b = {}, r_2 = {}".format(a,b,cur_r_2))
	print("Количество сопоставленных данных: {}".format(len(energies_parrots)))

	if cur_r_2 > max_r_2 and len(energies_parrots) >= set_limit:
		max_r_2 = cur_r_2
		best_shift = shift
		num_matched = len(energies_parrots)

	return(best_shift, max_r_2, a, b, num_matched)

def fit_and_check_shift(fit_function, energies_parrots, maxima_new, best_shift, shift, max_r_2, num_matched, set_limit, bounds = (np.array([1e-5, -np.inf]), np.array([np.inf, np.inf]))):

	a, b, cur_r_2 = sps.linregress(energies_parrots, maxima_new)[0:3]
	print("Fit and check function: a = {}, b = {}, r_2 = {}".format(a,b,cur_r_2))
	print("Количество сопоставленных данных: {}".format(len(energies_parrots)))

	if cur_r_2 > max_r_2 and len(energies_parrots) >= set_limit:
		max_r_2 = cur_r_2
		best_shift = shift
		num_matched = len(energies_parrots)

	return(best_shift, max_r_2, a, b, num_matched)

def check_best_shift(cur_r_2, max_r_2, energies_parrots, set_limit, shift, best_shift):

	if cur_r_2 > max_r_2 and len(energies_parrots) >= set_limit:
		max_r_2 = cur_r_2
		best_shift = shift
		num_matched = len(energies_parrots)

	return(best_shift, max_r_2, num_matched)
"""

def calc_av_squared_deviation_from_trend(fit_function, a, b, X, Y):
	'''
	Calculate averaged squared deviation of a data set from linear trend (sigma). Only y deviation is considered.

	fit_function --- fitting function (linear).
	a, b --- slope and intercept
	X, Y --- arrays of X and Y data values.
	'''

	dev = np.std(Y - fit_function(X, a, b))

	return(dev)

def f_rem_out(fit_function, foldername_ac, a, b, a1, b1, energies_parrots, maxima_new, is_zero, N_out=5.0):
	av_dev_lim = N_out*calc_av_squared_deviation_from_trend(fit_function, a, b, energies_parrots, maxima_new)
	energies_parrots_rem_out = np.zeros_like(energies_parrots)
	maxima_new_rem_out = np.zeros_like(maxima_new)
	for i in range(len(energies_parrots)):
		if abs(maxima_new[i] - fit_function(energies_parrots[i], a, b)) < av_dev_lim:
			energies_parrots_rem_out[i] = energies_parrots[i]
			maxima_new_rem_out[i] = maxima_new[i]
	energies_parrots_rem_out = energies_parrots_rem_out[np.abs(energies_parrots_rem_out) > is_zero]
	maxima_new_rem_out = maxima_new_rem_out[np.abs(maxima_new_rem_out) > is_zero]

	#Построение нового графика (без выпадающих точек).
	calibration_file = foldername_ac+"/Калибровка_same_rem_out" + ".png"
	dp.en_wf_plot(energies_parrots_rem_out, maxima_new_rem_out, calibration_file, style = '.')
	##########NEW##############
	plt.figure(figsize=(10.5, 9.0), dpi=200)
	plt.plot(energies_parrots_rem_out, maxima_new_rem_out, '.', color='k', lw=1.5)
	plt.plot([np.amin(energies_parrots_rem_out),np.amax(energies_parrots_rem_out)], [a*np.amin(energies_parrots_rem_out)+b, a*np.amax(energies_parrots_rem_out)+b], '-', lw=1.5, label = '1')
	plt.plot([np.amin(energies_parrots_rem_out),np.amax(energies_parrots_rem_out)], [a1*np.amin(energies_parrots_rem_out)+b1, a1*np.amax(energies_parrots_rem_out)+b1],  '-', lw=1.5, label = '2')
	plt.grid()
	plt.legend()
	plt.savefig(calibration_file, bbox_inches='tight')
	plt.close()

def parse_cmd_line():
	parser = argparse.ArgumentParser(description='Match the data with corresponding pulse energies.')
	parser.add_argument('folder_en', type=str, help="path to the folder containing energy files")
	parser.add_argument('folder_ac', type=str, help="path to the folder containing subfolders with data to be matched")
	parser.add_argument('ft', type=str, help="type of the data to be matched. Possible values: ac, mode, spectr, interf")
	parser.add_argument('ext', type=str, help="file extension. Possible values: .dat, .png, .tif, .bin")
	parser.add_argument('--old_osc', help="read data as for old version of oscilloscope (ms count). Otherwise, read as for new version (use time stamp with days and hours).", action="store_true")
	parser.add_argument('--try', dest='try_no_losses', help="check if there were no lost files, and if so, compare them one by one (without trying to correct lost)", action="store_true")
	parser.add_argument('--same', dest='same_computer', help="for data obtained at the same computer. Comparation directly by time.", action="store_true")
	parser.add_argument('--inv', help="for inverted acoustics data", action="store_true")
	parser.add_argument('-ro','--rem_out', dest='rem_out', help="remove points that outlines trend more than N_out sigma values.", action="store_true")
	parser.add_argument('--ac_lims', type=float, nargs=2, help="remove points that outlines trend more than N_out sigma values.")

	args = parser.parse_args()
	args.ac_lims = np.asarray(args.ac_lims)

	return(args)

def main():
	#%% Константы
	filetype_en = ".dat"
	shift_min = -3
	shift_max = 5
	set_limit = 10 # Минимальная длина выборки, начиная с которой данные считаются успешно сопоставленными.
	use_fon = True # Использовать ли сигнал "фона" (земли) с диода.
	col_fon = 6 #"старая" версия - 8, "новая" - 6.
	col_trig = 8 #"старая" версия - 6, "новая" - 8.
	col_en = 9 # номер столбца с энергиями.
	use_trig = True # Использовать ли колонку с сигналом строба, или считать начало записи выборки в начале файла с энергиями.
	r2_coeff = 1.0 # Если в режиме без коррекции пропусков (--try) r^2 > r2_coeff*r^2 в обычном режиме сопоставления, то используется режим без коррекции пропусков.
	low_r2_threshold = 0.1 # Если r^2 ниже этого значения, сопоставление считается проведённым неуспешно, выполняется попытка автоматического поиска shift (по начальным временам сопоставляемых данных), и сопоставление при этом значении shift.
	bounds_inv = (np.array([-np.inf, -np.inf]), np.array([-1e-5, np.inf])) #Границы подбора параметров аппроксимации для "перевёрнутой" акустики.
	bounds_not_inv = (np.array([1e-5, -np.inf]), np.array([np.inf, np.inf])) #Границы подбора параметров аппроксимации для "не перевёрнутой" акустики.
	DELTA = 10 #[мс] Максимально допустимый сдвиг текущего кадра по времени относительно рассчётного времени. По умолчанию DELTA = 15 мс.
	is_zero = 1e-8 #Условие "обращения в ноль" значений float при проверках.
	is_zero_int = 1e-1 #Условие "обращения в ноль" значений float при проверках, когда сравниваются целые числа.
	N_out = 5.0 #Если данные отклоняются от линейной зависимости более, чем на N_out*sigma, то они выбрасываются из завершающего сопоставления.

	#%% Начальные параметры
	#no_gaps = False #Используется обычный метод сопоставления (с коррекцией пропусков) - False, или без коррекции пропусков (в режиме --try). Флаг, используемый в программе для вывода в файл режима сопоставления, начальное значение менять не нужно.

	print("\n!!!")
	print("use_fon = {}".format(use_fon))
	print("use_trig = {}".format(use_trig))
	print(f'DELTA = {DELTA}\n')

	# Parse command line arguments.
	args = parse_cmd_line()
	print(args)

	folder_en = args.folder_en
	folder_ac = args.folder_ac
	ft = args.ft
	ext = args.ext

	print(folder_en)

	print("Use try_no_losses={}, same_computer={}, old_osc={}, inv={}".format(args.try_no_losses, args.same_computer, args.old_osc, args.inv))

	f_params = open(os.path.join(folder_ac, "Calibration_parameters.txt"), 'w')
	f_params.write("foldername_ac \t best_shift \t max_r^2 \t Number of the matched files\n")

	if ext=='.dat': #Задание области, в которой будет производиться интегрирование (для люминесценции)
		area = list(int(element) for element in input("Введите границы области интегрирования, в которой наблюдается сигнал. По умолчанию: 0 1920 0 1200\n").split())
	else:
		area=(0,1920,0,1200) #По умолчанию.

	if len(area) == 0:
		area = (0,1920,0,1200)
		print(f"Использую значение по умолчанию. area = {area}")
	else:
		print(f"area = {area}")

	#Путь к папке для копирования сопоставленных мод.
	folder_to_write = os.path.abspath(os.path.join(folder_ac, os.pardir))+"/"+folder_ac.split("/")[-2] + "_matched" + os.sep

	#Удаление папки, если она существует, и создание её заново.
	if os.path.exists(folder_to_write):
		sh.rmtree(folder_to_write)
	os.makedirs(folder_to_write)

	#Список файлов с энергиями и файлов с акустикой, сортировано в алфавитном порядке.
	filenames_en = [os.path.join(folder_en, f) for f in sorted(os.listdir(folder_en)) if f.endswith(filetype_en)]
	foldernames_ac = sorted(next(os.walk(folder_ac))[1])
	#print(filenames_en)
	#print(foldernames_ac)

	for filename_en in filenames_en:
		#Название файла с энергиями (без расширения и пути).
		filename_en_splitted = filename_en.split(os.sep)[-1]
		filename_en_splitted = filename_en_splitted.split(filetype_en)[0]

		for foldername_ac in foldernames_ac:
			foldername_ac_splitted = foldername_ac.split(os.sep)[-1]
			foldername_ac = os.path.join(folder_ac, foldername_ac)

			if filename_en_splitted == foldername_ac_splitted:
				print(foldername_ac)
				time_en, energies, i_start, lc, signal, fon = dp.read_en_all_data(filename_en, col_en=col_en, col_fon=col_fon, col_trig=col_trig)
				graph_file = foldername_ac+"/Фон.png"
				fon_i = np.linspace(0, len(fon), len(fon))
				dp.en_wf_plot(fon_i, fon, graph_file, style = '-')

				#i_start = 0
				if args.old_osc:
					filenames_ac_times, filenames_ac_info = dp.make_file_list_to_compare(foldername_ac, ext)
				else:
					filenames_ac_times, filenames_ac_info = dp.make_file_list_to_compare_new_program(foldername_ac, ext)

				if len(filenames_ac_times) == 0:
					print("В папке с сопоставляемыми данными нет файлов искомого типа.")
					break

				###TEMP####
				print("\n!!!")
				print(foldername_ac_splitted)
				###TEMP####
				if 'kalibr' in foldername_ac_splitted or 'calibr' in foldername_ac_splitted:
					maxima = dp.read_maxima(foldername_ac, filenames_ac_times, filenames_ac_info, ext, area=area, fon_coeff=1.0, old_osc=args.old_osc, limit_max = True, inv=False, ac_lims=None)
					print("HERE1!!!\n")
				else:
					maxima = dp.read_maxima(foldername_ac, filenames_ac_times, filenames_ac_info, ext, area=area, fon_coeff=1.0, old_osc=args.old_osc, limit_max = False, inv=args.inv, ac_lims=args.ac_lims)
					print("HERE2!!!\n")

				if ft == 'ac' and args.same_computer:
					calibration = dp.compare_not_save_same_computer(filename_en, foldername_ac, ext=ext, col_en=col_en, col_fon=col_fon, col_trig=col_trig, use_fon=use_fon, DELTA=DELTA)
					### Подготовка данных для графика.
					j_start = 0
					energies_parrots = []
					maxima_new = []
					for i in range(0, len(filenames_ac_times)):
						for j in range(j_start, calibration.shape[0]):
							if abs(filenames_ac_times[i] - calibration[j, 1]) < is_zero_int and abs(calibration[j,0]) > is_zero and abs(maxima[i]) > is_zero:
								energies_parrots.append(calibration[j,0])
								maxima_new.append(maxima[i])
								j_start = j
								break

					energies_parrots = np.array(energies_parrots)
					maxima_new = np.array(maxima_new)

					### Построение графика (для случая "same_computer").
					calibration_file = foldername_ac+"/Калибровка_same" + ".png"
					dp.en_wf_plot(energies_parrots, maxima_new, calibration_file, style = '.')

					#########NEW#############
					a, b, cur_r_2 = sps.linregress(energies_parrots, maxima_new)[0:3]
					a1, b1 = sps.linregress(energies_parrots, maxima_new)[0:2]
					f_params.write(foldername_ac+" \t --- \t {:.4f} \t {:.0f} \t 'Same' method.\n".format(cur_r_2, len(energies_parrots)))

					if args.rem_out:
						f_rem_out(fit_function, foldername_ac, a, b, a1, b1, energies_parrots, maxima_new, is_zero, N_out=N_out)
						###########END_NEW#########
					#!!! Пока только строится новый график! Реально данные не исправляются.
					#########################

					#Сопоставление (запись в файл).
					folder_ac_to_write = folder_to_write + filename_en_splitted + os.sep
					os.makedirs(folder_ac_to_write)
					dp.compare_same_computer(filename_en, foldername_ac, folder_ac_to_write, ext=ext, col_en=col_en, col_fon=col_fon, col_trig=col_trig, use_fon=use_fon, DELTA=DELTA)
				else:
					best_shift = -1
					max_r_2 = -1
					num_matched = 0
					a = 0; b = 0
					best_shift_old = -1
					max_r_2_old = -1
					num_matched_old = 0
					a_old = 0; b_old = 0

					for shift in range(shift_min, shift_max+1):
						energies_parrots = []
						maxima_new = []
						print("shift = {}".format(shift))

						if args.try_no_losses == True and dp.check_if_there_were_lost(filename_en, foldername_ac, ext, col_en=col_en, col_fon=col_fon, col_trig=col_trig, use_fon = use_fon, old_osc=args.old_osc):
							energies_parrots_old = []
							maxima_new_old = []
							calibration = dp.compare_not_save_no_lost(filename_en, foldername_ac, ext, shift = shift, col_en=col_en, col_fon=col_fon, col_trig=col_trig, use_fon = use_fon, use_trig=use_trig, old_osc=args.old_osc, DELTA=DELTA) # No gaps
							# old -> general method (with gaps, old)
							calibration_old = dp.compare_not_save_new_method(filename_en, foldername_ac, ext, shift = shift, col_en=col_en, col_fon=col_fon, col_trig=col_trig, use_fon = use_fon, use_trig=use_trig, old_osc=args.old_osc, DELTA=DELTA) #Для согласованности надо не забыть поменять метод ниже.
							#For general (old) method.
							j_start = 0
							for i in range(0, len(filenames_ac_times)):
								for j in range(j_start, calibration_old.shape[0]):
									if (abs(filenames_ac_times[i] - calibration_old[j, 1]) < is_zero_int) and (abs(calibration_old[j,0]) > is_zero) and (abs(maxima[i]) > is_zero):
										energies_parrots_old.append(calibration_old[j,0])
										maxima_new_old.append(maxima[i])
										j_start = j
										break
							energies_parrots_old = np.array(energies_parrots_old)
							maxima_new_old = np.array(maxima_new_old)
							if energies_parrots_old.size <= 1:
								print('Size of the matched energies array <= 1 for shift {}'.format(shift))
								continue

							a_old, b_old, cur_r_2_old = sps.linregress(energies_parrots_old, maxima_new_old)[0:3]
							if cur_r_2_old > max_r_2_old and len(energies_parrots_old) >= set_limit:
								max_r_2_old = cur_r_2_old
								best_shift_old = shift
								num_matched_old = len(energies_parrots_old)
						else:
							calibration = dp.compare_not_save_new_method(filename_en, foldername_ac, ext, shift = shift, col_en=col_en, col_fon=col_fon, col_trig=col_trig, use_fon = use_fon, use_trig=use_trig, old_osc=args.old_osc, DELTA=DELTA) #Для согласованности надо не забыть поменять метод ниже.

						# Калибровка для нового метода (если используется флаг --try и не были пропущены стробы.
						# Либо калибровка для старого (general) метода, если используется только он.
						j_start = 0
						for i in range(0, len(filenames_ac_times)):
							for j in range(j_start, calibration.shape[0]):
								if abs(filenames_ac_times[i] - calibration[j, 1]) < is_zero_int and abs(calibration[j,0]) > is_zero and abs(maxima[i]) > is_zero:
									energies_parrots.append(calibration[j,0])
									maxima_new.append(maxima[i])
									j_start = j
									break

						energies_parrots = np.array(energies_parrots)
						maxima_new = np.array(maxima_new)

						#%% Построение графиков.
						if args.try_no_losses == True and dp.check_if_there_were_lost(filename_en, foldername_ac, ext, col_en=col_en, col_fon=col_fon, col_trig=col_trig, use_fon = use_fon, old_osc=args.old_osc):
							# No gaps and general method.
							calibration_file = foldername_ac+"/Калибровка_" + str(shift) + "_no_gaps.png"
							dp.en_wf_plot(energies_parrots, maxima_new, calibration_file, style = '.', dpi=200)
							calibration_file = foldername_ac+"/Калибровка_" + str(shift) + ".png"
							dp.en_wf_plot(energies_parrots_old, maxima_new_old, calibration_file, style = '.', dpi=200)
						else:
							# General method.
							calibration_file = foldername_ac+"/Калибровка_" + str(shift) +".png"
							dp.en_wf_plot(energies_parrots, maxima_new, calibration_file, style = '.', dpi=200)

						if energies_parrots.size <= 1:
							print('Size of the matched energies array <= 1 for shift {}'.format(shift))
							continue

						a, b, cur_r_2 = sps.linregress(energies_parrots, maxima_new)[0:3]
						if cur_r_2 > max_r_2 and len(energies_parrots) >= set_limit:
							max_r_2 = cur_r_2
							best_shift = shift
							num_matched = len(energies_parrots)
						best_shift, max_r_2, a, b, num_matched = fit_and_check_shift(fit_function, energies_parrots, maxima_new, best_shift, shift, max_r_2, num_matched_old, set_limit, bounds=bounds_inv)

					#%% Если максимальный r2 меньше порога, запускаем автоматический поиск "правильного" сдвига.
					if args.try_no_losses == True and dp.check_if_there_were_lost(filename_en, foldername_ac, ext, col_en=col_en, col_fon=col_fon, col_trig=col_trig, use_fon = use_fon, old_osc=args.old_osc):
						if max_r_2 < low_r2_threshold and max_r_2_old < low_r2_threshold:
							best_shift = dp.shift_search(filenames_ac_times, dt = 0.1, shift_border_min = -3, shift_border_max = 5)
							shift_corr = True
							print("Best_shift was recalculated automatically.")
							print("New best_shift = {}".format(best_shift))
						else:
							shift_corr = False
					else:
						if max_r_2 < low_r2_threshold:
							best_shift = dp.shift_search(filenames_ac_times, dt = 0.1, shift_border_min = -3, shift_border_max = 5)
							shift_corr = True
							print("Best_shift was recalculated automatically.")
							print("New best_shift = {}".format(best_shift))
						else:
							shift_corr = False

					# осталось запустить разметку интерферограмм с нужным best_shift
					folder_ac_to_write = os.path.join(folder_to_write, filename_en_splitted)
					os.makedirs(folder_ac_to_write)
					if args.try_no_losses == True and dp.check_if_there_were_lost(filename_en, foldername_ac, ext, col_en=col_en, col_fon=col_fon, col_trig=col_trig, use_fon = use_fon, old_osc=args.old_osc):
						print('max r_2 = {}, max r_2 old = {}, optimal_shift = {}, optimal_shift_old = {}'.format(max_r_2, max_r_2_old, best_shift, best_shift_old))
						if max_r_2 > max_r_2_old * r2_coeff:
							print("Используется метод без пропусков.")
							dp.compare_no_lost(filename_en, foldername_ac, folder_ac_to_write, ext=ext, shift = best_shift, use_fon = use_fon, col_en=col_en, col_fon=col_fon, col_trig=col_trig, use_trig=use_trig, old_osc=args.old_osc, DELTA = DELTA)
							no_gaps = True
						else:
							print("Используется обычный метод сопоставления.")
							dp.compare_new_method(filename_en, foldername_ac, folder_ac_to_write, ext=ext, shift = best_shift_old, use_fon = use_fon, col_en=col_en, col_fon=col_fon, col_trig=col_trig, use_trig=use_trig, old_osc=args.old_osc, DELTA = DELTA) #Для согласованности надо не забыть поменять метод выше.
							no_gaps = False
					else:
						print('max r_2: {}, optimal_shift = {}'.format(max_r_2, best_shift))
						dp.compare_new_method(filename_en, foldername_ac, folder_ac_to_write, ext=ext, shift = best_shift, use_fon = use_fon, col_en=col_en, col_fon=col_fon, col_trig=col_trig, use_trig=use_trig, old_osc=args.old_osc, DELTA = DELTA) #Для согласованности надо не забыть поменять метод выше.
						no_gaps=False

					if shift_corr:
						f_params.write(foldername_ac+" \t {:.0f} \t {:.4f} \t {:.0f} \t Low r^2. Shift was found automatically\n".format(best_shift, max_r_2, num_matched))
					elif no_gaps:
						f_params.write(foldername_ac+" \t {:.0f} \t {:.4f} \t {:.0f} \t no gaps method\n".format(best_shift, max_r_2, num_matched))
					else:
						f_params.write(foldername_ac+" \t {:.0f} \t {:.4f} \t {:.0f} \t general method\n".format(best_shift, max_r_2, num_matched))

				#calibr_file = os.path.join(folder_ac, foldername_ac, "Calibr_" + str(shift) + ".png")
				#dp.en_wf_plot(energies_parrots, maxima_new, calibr_file, style = '.')

	f_params.close()

if __name__=='__main__':
	main()
