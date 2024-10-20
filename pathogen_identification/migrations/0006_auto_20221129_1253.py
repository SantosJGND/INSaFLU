# -*- coding: utf-8 -*-
# Generated by Django 1.11.29 on 2022-11-29 12:53
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('pathogen_identification', '0005_auto_20221128_0915'),
    ]

    operations = [
        migrations.AddField(
            model_name='rundetail',
            name='depleted_reads',
            field=models.IntegerField(blank=True, default=0, null=True),
        ),
        migrations.AddField(
            model_name='rundetail',
            name='depleted_reads_percent',
            field=models.FloatField(blank=True, default=0, null=True),
        ),
        migrations.AddField(
            model_name='rundetail',
            name='enriched_reads',
            field=models.IntegerField(blank=True, default=0, null=True),
        ),
        migrations.AddField(
            model_name='rundetail',
            name='enriched_reads_percent',
            field=models.FloatField(blank=True, default=0, null=True),
        ),
    ]
