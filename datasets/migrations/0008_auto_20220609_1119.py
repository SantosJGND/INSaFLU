# -*- coding: utf-8 -*-
# Generated by Django 1.11.18 on 2022-06-09 11:19
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('datasets', '0007_auto_20220608_1445'),
    ]

    operations = [
        migrations.AddField(
            model_name='datasetconsensus',
            name='name',
            field=models.CharField(blank=True, db_index=True, max_length=200, null=True, verbose_name='Name'),
        ),
        migrations.AddField(
            model_name='datasetconsensus',
            name='type_subtype',
            field=models.CharField(blank=True, max_length=50, null=True),
        ),
    ]