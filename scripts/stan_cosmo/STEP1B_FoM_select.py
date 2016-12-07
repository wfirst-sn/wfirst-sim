from numpy import *
import matplotlib.pyplot as plt

def read_csv(fl):
    f = open(fl)
    lines = f.read().split('\n')
    f.close()

    headers = lines[0].split(",")

    data = {}
    for item in headers:
        data[item] = []
    
    for line in lines[1:]:
        parsed = line.split(",")
        if len(parsed) > 1:
            for i in range(len(headers)):
                try:
                    data[headers[i]].append(float(parsed[i]))
                except:
                    data[headers[i]].append(parsed[i])


    data["maxz"] = []
    for i in range(len(data["0.4<z<0.8"])):
        if data["1.7<z<2.0"][i] > 0:
            data["maxz"].append(2.0)
        elif data["1.4<z<1.7"][i] > 0:
            data["maxz"].append(1.7)
        elif data["1.1<z<1.4"][i] > 0:
            data["maxz"].append(1.4)
        elif data["0.8<z<1.1"][i] > 0:
            data["maxz"].append(1.1)
        elif data["0.4<z<0.8"][i] > 0:
            data["maxz"].append(0.8)

    data["Spectra S/N"] = [float(item.split("+")[-1]) for item in data["Spectra S/N"]]
    data["Has Ground"] = [float(eval(item)) for item in data["Has Ground"]]

    return data


def eval_inds(inds, data, keys):
    metric = 0
    for key in keys:
        vals = unique(data[key])
        frac = 1./len(vals)

        for val in vals:

            observed_frac = float(sum(data[key][inds] == val))/float(len(data[key][inds]))
            metric -= log(observed_frac/frac)**2.
    return metric


#survey_00004,6:160720.txt,0.050,206.140,shallow+medium+deep+shallow+shallow,13.50+23.20+38.60+45.00,True,28.8,0x0,5.0,300x2,676,2464,340,322,0,0,

data = read_csv("summary.csv")

nsurveys = 60

print data.keys()

print len(data["Time Used"])

print median(data["Time Used"])

print data["maxz"]

for i in range(len(data["Time Used"]))[::-1]:
    bad_survey = 0

    if float(data["Time Used"][i]) > 230 or float(data["Time Used"][i]) < 205 or float(data["Spectra S/N"][i]) < 20:
        bad_survey = 1

    
    """
    print "HACK!!!!!!!"
    if data["Spectra S/N"][i] == 40:
        bad_survey = 1

    """

    """
    print "HACK!!!!!!!"
    if data["Wide Square Degrees"][i] <40:
        bad_survey = 1
        """


    if bad_survey:
        for key in data.keys():
            del data[key][i]
    



    

print len(data["Time Used"])

for key in data.keys():
    data[key] = array(data[key])


best_metric = -1000000

for attempt in range(200000):
    if attempt % 5000 == 0:
        print attempt
        #print "HACK!!!!!"
    inds = random.choice(arange(len(data["Time Used"])), size = nsurveys, replace = False) # HACK!!!!!!!!!!!!!

    #print inds
    assert len(unique(inds)) == len(inds)

    keys = ["Deep Square Degrees", "Spectra S/N", "Has Ground", "Pixel Scale", "maxz"]

    metric = eval_inds(inds, data, keys)
    if metric > best_metric:
        print "Better ", metric
        print "attempt ", attempt
        best_metric = metric
        best_inds = inds

inds = best_inds

print "survey_tests_new/" + "  survey_tests_new/".join(list(data["Survey"][inds]))

print "survey_tests_limited/" + "  survey_tests_limited/".join(list(data["Survey"][inds]))


plt.figure(figsize = (15,15))
for i in range(len(keys)):
    for j in range(len(keys)):
        plt.subplot(len(keys), len(keys), 1+i+j*len(keys))
        if i != j:
            plt.plot(data[keys[i]][inds], data[keys[j]][inds], '.', color = 'b')
        else:
            plt.hist(data[keys[i]][inds])
        plt.title(keys[j] + " vs " + keys[i])
plt.savefig("scatter.pdf")
plt.close()
