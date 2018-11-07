# -*- coding: utf-8 -*-
__author__ = 'qsong'


def test_var_args(arg1, arg2, arg3, *args):
    print "arg1", arg1
    print "arg2", arg2
    print "arg3", arg3
    for arg in args:
        print "other args:", arg


def test_var_kwargs(farg, **kwargs):
    print "formal arg:", farg
    for key in kwargs:
        print "another keyword arg: %s: %s" % (key, kwargs[key])


if __name__ == "__main__":

    arg_dict = {
        "myarg2": 4,
        "suibian": 7
    }

    test_var_args(1, 2, 3, "a", True, "HHH")
    test_var_kwargs(2, myarg2="two", myarg3=3)
    test_var_kwargs(2, arg_dict)


