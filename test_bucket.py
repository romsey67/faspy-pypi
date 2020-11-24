from faspy.basel.liquidity.buckets import bucketing
from faspy.interestrate import discount_curve as dcurve

# %%

rate = {"value_date": "2020-10-30", "st_busday": "Modified Following",
        "st_ratebasis": "Money Market", "st_daycount": "Actual/365",
        "lt_busday": "No Adjustment", "lt_frequency": "Semi-Annual",
        "lt_daycount": "Actual/Actual",
        "rates": {'O/N': 2.30, '1W': 2.35, '1M': 2.45, '3M': 2.55,
                  '6M': 2.65, '12M': 2.75, '1Y': 2.70, '2Y': 2.80, '3Y': 2.90,
                  '5Y': 3.00, '10Y': 3.10, '30Y': 3.25}}

df = dcurve.discount_factor_gen(rate, return_type="time")
data = [{"id": 0, "date": "2021-10-01", "amount": 100},
{"id": 1, "date": "2022-05-01", "amount": 290_000}]


bucket = [
    {"name": ">10Y and <=30Y", "from": 10, "to": 30},
    {"name": "<=1W", "from": 0, "to": 7/365},
    {"name": ">1W and <=1M", "from": 7/265, "to": 1/12},
    {"name": ">1M and <=3M", "from": 1/12, "to": 3/12},
    {"name": ">3M and <=6M", "from": 3/12, "to": 6/12},
    {"name": ">6M and <=1Y", "from": 6/12, "to": 1},
    {"name": ">1Y and <=2Y", "from": 1, "to": 2},
    {"name": ">2Y and <=3Y", "from": 2, "to": 3},
    {"name": ">3Y and <=5Y", "from": 3, "to": 7},
    {"name": ">5Y and <=10Y", "from": 5, "to": 10},
    {"name": " >30Y", "from": 30, "to": 1_000}
          ]

newdata = bucketing(data, bucket, "2020-01-01", dfcurve=df)
print(newdata)
