import datetime


class TimeTool:
    @staticmethod
    def date2mjd(date):
        t0 = datetime.date(1858, 11, 17)
        mjd = (date - t0).days
        return mjd

    @staticmethod
    def date2yyyymmdd(date):
        year_str = str(date.year)
        month_str = str(date.month).rjust(2, '0')
        day_str = str(date.day).rjust(2, '0')
        return year_str + month_str + day_str

    @staticmethod
    def mjd2yyyymmdd(mjd):
        date = TimeTool.mjd2date(mjd)
        return TimeTool.date2yyyymmdd(date)

    @staticmethod
    def mjd2date(mjd):
        t0 = datetime.date(1858, 11, 17)
        return t0 + datetime.timedelta(days=mjd)

    @staticmethod
    def isLeap(year):
        return True if (year % 4 == 0 and year % 100 != 0) else False

    @staticmethod
    def year_frac(date: datetime.date):
        year = date.year
        if TimeTool.isLeap(year):
            days = [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366]
        else:
            days = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]

        month = date.month
        day = date.day

        return year + (days[month - 1] + day) / days[12]

    @staticmethod
    def yyyymm(date: datetime.date):
        year = str(date.year)
        month = str(date.month).rjust(2, '0')
        return year + month

    @staticmethod
    def list_of_unsed_days(begin, end, unused_days: list):
        used_days = []
        max_days = (end - begin).days
        for i in range(max_days + 1):
            this_date = begin + datetime.timedelta(days=i)
            if this_date not in unused_days:
                used_days.append(this_date)

        return used_days

    @staticmethod
    def average_date_of_series(series: list):
        mjd = 0
        for iter in range(len(series)):
            mjd += TimeTool.date2mjd(series[iter])
        return TimeTool.mjd2date(mjd / len(series))
