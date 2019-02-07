есть два варианта запуска, в обоих нужен файл в формате фаста (пики, файл *.seq), имя файла = аргумент программы
также всегда являются аргументами 
а) максимальная длина спейсера (по умолчанию я ставлю всегда 30, внутри программы ограничение до 100)
б) путь к папке где лежит полногеномный файл (для человека мыши ups2kb.plain, для арабидопсиса ups1500.plain, они лежат в папках 
matrices\hocomoco_hs402,  matrices\hocomoco_mm358,  matrices\dapseq528 = для человека, мыши и арабидопсиса
(также в папках лежат партнерские мотивы в формате MEME, они будут нужны когда будем прикручивать Tomtom, пока не нужны)

вариант запуска 1: anchor_vs_many_partners

партнерские мотивы задаются библитеками 
а) Hocomoco (человек или мышь, CORE или FULL collection) отсюда http://hocomoco11.autosome.ru/downloads_v11
б) Dapseq (арабидопсис) отсюда http://neomorph.salk.edu/dev/pages/shhuang/dap_web/pages/browse_table_aj.php

якорный мотив в формате *.motif (матрица частот) обязан быть в папке где программа будет считать, имя якорного мотива задано hocomoco0.motif или dapseq0.motif
для баз Hocomoco или DApseq


вариант запуска 2: anchor_vs_one_partner

якорный и партнерский мотив = аргументы программы

командные строки запуска - в папке run
самое нужное из выходных файлов - два файла out_hist, out_pval - в папке out


файл out_hist = сводная таблица значимостей (далее объяснение для варианта запуска anchor_vs_many_partners, для anchor_vs_one_partner в целом аналогично)

в колонках:


MotifNumber	MotifName - номер и название партнерского мотива

Значимости композиционного элемента
AnyPv	FullPv	PartPv	OverPv	SpacPv	:	Значимости p-value сопредставленности мотивов : (перекрывание или спейсер), полное перекрывание, частичное перекрывание, (полное или частичное перекрывание), спейсер


Асимметрия консервантивности мотивов в составе композиционного элемента

в программе используется по пять порогов любой весовой матрицы, поэтому рассчитвается асимметрия в матрице Log10(p-value) по матрице 5х5 порог_якоря VS порог партнера
AnyAssy	FullAssy	PartAssy	OverAssy	SpacAssy  :  коэффициенты асимметрии


Сравнение мотивов

Motif_comparison производит сравнение мвтрицы якоря и партнера и результат (p-value) записывается в колонку SimTot, значение < 0.05 доказывает что матрицы значимо похожи

далее записываются колонки с предварительными расчетами p-value по Евклидовому и Пирсоновскому расстоянию межлу матрицами (с малым числом итераций), SimEDpre	SimPCpre	
если эти тесты дают значимость - тогда считаются тесты с большим числом итераций, SimED	SimPC
