# -*- coding: utf-8 -*-
# Generated by Django 1.11.18 on 2020-05-04 19:33
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('settings', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='parameter',
            name='not_set_value',
            field=models.CharField(blank=True, db_index=True, max_length=50, null=True),
        ),
        migrations.AlterField(
            model_name='parameter',
            name='description',
            field=models.CharField(default='', max_length=500),
        ),
    ]
